// This stan model aims at inferring transmission strength between admin units
// based on cholera incidence and phylogenetic groups

// Forward/backward algorithm adapted from https://github.com/luisdamiano/stancon18/blob/master/stan/hmm_gaussian.stan
// and code generously provided from Claire (Justin's JHU PhD student)

functions {
  
  real prob_v0_calc(int n, real epsilon, int V, row_vector probs){
    real pplus = 1;
    
    for(v in 1:V){
      pplus *= 1-probs[v];
    }

    if(pplus<epsilon){
      pplus = epsilon;
    }else if(pplus==1){
      pplus = 1-epsilon;
    }
    return pplus;
  }
  
}

data {
  //Indexing
  int<lower=1> M;    // number of spatial locations (countries)
  int<lower=1> T;    // maximum number of time periods (years)
  int<lower=1> V;    // number of phylogenetic groups/strains
  int<lower=1> K;    // number of hidden states

  // Observed sequence data
  int<lower=1> N;    // number of space x time units 
  array[N] int<lower=0> N_sequences;    // total number of observed sequences across all strains (space x time)
  array[N,V] int<lower=0> y;    // number of observed sequences of each strain 
  
  // Observed Case Data
  array[N] int<lower=0> cases;    // number of cholera cases reported per location and time period
  
  // Population sizes
  vector[M] pop;
  
  // Mappings
  array[N] int<lower=0, upper=M> map_to_location;    // map from observation to location 
  array[N] int<lower=1, upper=T> map_to_time;        // map from observation to time period 
  array[M,T] int map_to_n;    //map from location/time back to N

  // Introductions
  matrix[M,V] intro_times0;    //upper bound on introduction times
  int<lower=1> N_gamma;    // number of strains with introductions
  int<lower=1> N_loc_gamma;    // number of locations with introductions
  int<lower=1> N_yr_gamma;    // number of years with introductions
  array[V] int map_v_to_gamma;
  array[M, V] int map_gamma_loc;
  array[T, V] int map_gamma_yr;
  
  real epsilon;    // small number to avoid log(0) issues
  
  // Spatial connectivity
  int<lower=0> N_edges;    // number of edges between locations
  array[M, M] int<lower=0> map_from_edge;    // map from edge number to origin/dest matrix
  array[N_edges] int <lower=1, upper=M> node1_by_dest;    // column 1 of the adjacency matrix
  vector[N_edges] dist_by_dest;    // column 1 of the adjacency matrix
  array[M] int<lower=0> start_by_dest;    // start of adjancency list per location
  array[M] int<lower=0> end_by_dest;      // end of adjancency list per location
  real p_w;
  vector[2] mu_w;    //spatial weights
  vector[2] sd_w;    //spatial weights

  // Priors
  vector[2] mu_tau; // population parameter - gravity
  real sd_tau; // population parameter - gravity
  real mu_zeta; // distance parameter - gravity
  real sd_zeta; // distance parameter - gravity
  real mu_kappa; // gravity constant
  real sd_kappa; // gravity constant
  real mu_delta; // persistence constant
  real sd_delta; // persistence constant
  real sd_lambda;    //strain-specific sequence rate
  real mu_eta; // case transformation
  real sd_eta; // case transformation
  real mu_e; // case under-reporting
  real sd_e; // case under-reporting
  array[N_gamma] vector[N_loc_gamma] p_gamma_loc;
  array[N_gamma] vector[N_yr_gamma] p_gamma_yr;
  real l_gamma;

}


parameters {
  //connectivity
  vector<lower=0.1, upper=0.6>[2] tau; // population parameter - gravity model
  real<lower = (1.5-mu_zeta)/sd_zeta, upper = (3-mu_zeta)/sd_zeta> zeta_std; // distance parameter - gravity model
  real<lower = 0.001, upper = 0.999> kappa; // log gravity constant
  vector<lower = -8, upper = 4> [N_edges] w;    // spatial connectivity parameter for each edge

  //persistence
  // vector<lower = -10, upper = 1>[M] log_delta;    // strain persistence from one year to the next in a given country
  real<lower = -10, upper = 1> log_delta;    // strain persistence from one year to the next in a given country

  //prevalence
  matrix<lower = log(1e-4), upper = log(max(N_sequences))>[N,V] log_lambda;    // for poisson approximation of prevalence

  //cases
  vector<lower = 0, upper = 1>[M] e_i;            // under-reporting of cholera cases (location-specific)
  real<lower = 0.2, upper = 0.8> eta; // case transformation
   
  //introductions
  array[N_gamma] vector<lower = 0, upper = 1>[N_loc_gamma] gamma_loc;
  array[N_gamma] vector<lower = 0, upper = 1>[N_yr_gamma] gamma_yr;
  
  //initial probabilities
  matrix<lower = 1e-9, upper = 1-1e-9>[M,V] pi_0;

}
 
transformed parameters {
  
  // connectivity
  real zeta = mu_zeta+sd_zeta*zeta_std; //distance
  vector[N_edges] xi; // connectivity parameter between different locations

  // vector[M] delta = exp(log_delta); //persistence
  real delta = exp(log_delta); //persistence
  
  // strain rate parameter & prevalence
  matrix[N,V] lambda = exp(log_lambda);    // sequencing rate parameter for poisson approximation of prevalence
  matrix[N,V] lambda_star;    // prevalence (based on sequencing)

  // introduction probability of strain v in country i at time t. This incorporates the decay function
  matrix[N,V] gamma;    //introduction from outside continent
  
  array[K] matrix<lower=0, upper=1>[N,V] p_col; // probability of colonization (transition)

  //Initialize matrices to track forward algorithm (alpha), backward algorithm (beta), & rho (alpha * beta)
  array[V] matrix[K,N] logalpha;    // log forward
  array[V] matrix<lower = 0, upper = 1>[K,N] alpha;    // forward
  
  {// connectivity
    for(i in 1:M){
      for (k in start_by_dest[i]:end_by_dest[i]) {// start other country (j) loop
          int j = node1_by_dest[k];    // origin location into m
          int ji = map_from_edge[j,i];
          xi[ji] = exp(log(kappa * (((pop[j])^tau[1]*(pop[i])^tau[2])/(dist_by_dest[ji]^zeta))) + w[ji]);
 
      }
    }
  }

  
  {// introductions
  for(v in 1:V){
    for(i in 1:M){
      for(t in 1:T){
        
        int n = map_to_n[i,t];
        real p_loc;
        real p_t;
        
        if(map_v_to_gamma[v]>0){
          if(map_gamma_loc[i,v]>0){
            p_loc = gamma_loc[map_v_to_gamma[v]][map_gamma_loc[i,v]];
          }else{
            p_loc = epsilon;
          }
          if(map_gamma_yr[t,v]>0){
            p_t = gamma_yr[map_v_to_gamma[v]][map_gamma_yr[t,v]];
          }else{
            p_t = epsilon;
          }
          gamma[n,v] = p_loc*p_t;
        }else{
          gamma[n,v] = 1e-5;
        }
        if(gamma[n,v]<1e-5){
          gamma[n,v] = 1e-5;
        }else if(gamma[n,v]>1-1e-5){
          gamma[n,v] = 1-1e-5;
        }
      }
    }
  }
  }// introductions
  
  
  {// start t = 1 probability of presence //
  
  // loop through strains
  for(v in 1:V){
    // loop through countries
    for(i in 1:M){
      vector[T] log_lambda_iv = log_lambda[map_to_n[i,1]:map_to_n[i,T],v]; // sequencing rate parameter for poisson approximation of prevalence
      array[T] int y_iv = y[map_to_n[i,],v]; // subset observations to country/strain
      int n_t1 = map_to_n[i,1]; // n at time 1

      for(k in 1:K){
        int z = k-1;
        real pi;
        if(z==1){
          pi = pi_0[i,v];
        }else{
          pi = 1-pi_0[i,v];
        }                    //probability at time 1 * probability of observation given state
        logalpha[v][k,n_t1] = log(pi) + log_obs_func(y_iv[1],log_lambda_iv[1],k-1,epsilon);
      }
      
      // convert 
      alpha[v][1:K,n_t1] = softmax(to_vector(logalpha[v][1:K,n_t1]));
    
      //initialize p_col for time 1 
      for(k in 1:K){
        p_col[k][n_t1,v] = epsilon;
      }

    
    }//end country loop
  }//end strain loop
  
  }// end t = 1 presence //
  
  {//prevalence at t = 1
  for(i in 1:M){
    int n_t1 = map_to_n[i,1]; // n at time 1
    real lambda_tot = 0;
      
    for(v in 1:V){
      lambda_tot += lambda[n_t1,v]*alpha[v][K,n_t1] + epsilon;
    }
      
    for(v in 1:V){
      lambda_star[n_t1,v] = (lambda[n_t1,v]*alpha[v][K,n_t1] + epsilon)/lambda_tot;
    }//end strain loop
  }//end country loop
    
  }//prevalence at t = 1
    

  {//// start FORWARD algorithm log p(z_t = j | x_{1:t})  ////
  
  // loop through time
  for(t in 2:T){  
    
    {// start colonization probability
      
    // loop over countries
    for(i in 1:M){
      // loop over strains
      for(v in 1:V){
        // mappings
        int n = map_to_n[i,t]; // n at t for country i
        int nm = map_to_n[i,t-1]; // n at t-1 for country i
        // initialize
        // real phi = 1-(1-exp(-gamma[nm,v])); // initialize rate of colonization & self-introduction
        real phi_base = 1-gamma[nm,v]; // initialize probability of colonization & self-introduction
        
        for (k in start_by_dest[i]:end_by_dest[i]) {// start other country (j) loop
          int j = node1_by_dest[k];    // origin location into m
          int nm_j = map_to_n[j,t-1]; // n at t-1 for country j
          int ji = map_from_edge[j,i];
          
          // update rate of colonization
          // prevalence*cases*connectivity ~ per case transmission for strain v multiplied by P(z_j=1)
          phi_base *= 1 - (1 - exp(-pow(lambda_star[nm_j, v] * cases[nm_j], eta) * xi[ji]));
        
        }// end other country loop
          
        // update rate of colonization to incorporate self re-introduction (delta)
        // previous state 0 -> current state 1
        p_col[1][n, v] = 1 - phi_base;

        // previous state 1 -> current state 1
        {
          real phi_present = phi_base;
          phi_present *= 1 - (1 - exp(-pow(lambda_star[nm, v] * cases[nm], eta) * delta));
          p_col[2][n, v] = 1 - phi_present;
        }
        // for(k in 1:K){
        //   int zprev = k-1;//state at time t-1
        //   if(zprev==1){
        //     // prevalence*cases*persistence ~ per case transmission for strain v
        //     phi *= 1-(1-exp(-((lambda_star[nm,v] * cases[nm])^eta * delta)));
        //     p_col[k][n,v] = 1-phi;//convert to pcol[2], this is transition from state 1->1, 1-pcol[2] is 1->0
        //   }else{
        //     p_col[k][n,v] = 1-phi;//convert to pcol[1], this is transition from state 0->1, 1-pcol[1] is 0->0
        //   }
        // }
      }// end strain loop
    }// end country loop
   
    }// end colonization probability
    
    {//transition probabilities & forward equation
    
    // loop over countries
    for(i in 1:M){
      // loop through strains
      for(v in 1:V){
        
        // initialize
        // matrix[K,K] A; // transition probabilities;
        // array[K] real accumulator; 
        // subset observed data to strain & country
        int y_iv = y[map_to_n[i,t],v];
        real log_lambda_iv = log_lambda[map_to_n[i,t],v]; // sequencing rate parameter for poisson approximation of prevalence
        // mappings
        int n = map_to_n[i,t]; // n at t for country i
        int nm = map_to_n[i,t-1]; // n at t-1 for country i
        
        // observation log probabilities
        real ll_absent;
        real ll_present;
        
        // transition probabilities
        real a01;  // 0 -> 1
        real a11;  // 1 -> 1
        real a00;  // 0 -> 0
        real a10;  // 1 -> 0
        
        // transition probabilities, clamped away from 0/1
        a01 = p_col[1][n, v];
        if (a01 == 1) a01 = 1 - epsilon;
        else if (a01 == 0) a01 = epsilon;

        a11 = p_col[2][n, v];
        if (a11 == 1) a11 = 1 - epsilon;
        else if (a11 == 0) a11 = epsilon;

        a00 = 1 - a01;
        if (a00 == 1) a00 = 1 - epsilon;
        else if (a00 == 0) a00 = epsilon;

        a10 = 1 - a11;
        if (a10 == 1) a10 = 1 - epsilon;
        else if (a10 == 0) a10 = epsilon;
        
        // observation log-probabilities
        if (y_iv == 0) {
          ll_absent = log1m(epsilon);
        } else {
          ll_absent = log(epsilon);
        }
        ll_present = poisson_log_lpmf(y_iv | log_lambda[n, v]);
      
        // observation process --> forward algorithm
        // state 0 at time t
        logalpha[v][1, n] = log_sum_exp(logalpha[v][1, nm] + log(a00), logalpha[v][2, nm] + log(a10)) + ll_absent;

        // state 1 at time t
        logalpha[v][2, n] = log_sum_exp(logalpha[v][1, nm] + log(a01), logalpha[v][2, nm] + log(a11) ) + ll_present;

        // convert 
        {
          real denom = log_sum_exp(logalpha[v][1, n], logalpha[v][2, n]);
          alpha[v][1, n] = exp(logalpha[v][1, n] - denom);
          alpha[v][2, n] = exp(logalpha[v][2, n] - denom);
        }
        
      }// end strain loop
    }// end country loop
    
    
    {//prevalence 
    for(i in 1:M){//country loop
      int n = map_to_n[i,t]; // n at t for country i
      real lambda_tot = 0;
      
      for(v in 1:V){//start strain loop
        lambda_tot += lambda[n,v]*alpha[v][2,n] + epsilon;
        
      }//end strain loop
      
      for(v in 1:V){
        
        lambda_star[n,v] = (lambda[n,v]*alpha[v][2,n] + epsilon)/lambda_tot;
        
      }//end strain loop
    }//end country loop
    
    }//prevalence 
    
    }// end transition probabilities & forward algorithm
  
  }//end time loop
  
  }// end FORWARD algorithm

     
}//end transformed parameters block


model {
  
  {// target alpha
    for(i in 1:M){
      for(v in 1:V){
          int n = map_to_n[i,T];
          target += log_sum_exp(logalpha[v][,n]);
        }
      }
  }// target alpha

  {//case & N_sequence obs
  matrix[N,V] alpha_present;
  for(n in 1:N){
    for(v in 1:V){
      alpha_present[n,v] = alpha[v][2,n];
    }
  }
  for(n in 1:N){
    real alpha_0 = prob_v0_calc(n,epsilon,V,alpha_present[n,]); //probability of any strains being present
    int i = map_to_location[n];

    //case observation
    if (cases[n] == 0) {
        target += log_sum_exp([log(alpha_0), log(1-alpha_0) + poisson_lpmf(0 | (1/e_i[i]) + N_sequences[n])]);
      } else {
        target += log(1-alpha_0) + poisson_lpmf(cases[n] | (cases[n]/e_i[i]));
      }
  }
  }//case_obs


  // --- Parameter priors ---

  {// case underreporting & sampling rate
  eta ~ beta(mu_eta*sd_eta, (1-mu_eta)*sd_eta);
  for(i in 1:M){
    //case underreporting
    e_i[i] ~ beta(mu_e*sd_e, (1-mu_e)*sd_e);

  }
  }//case underreporting & sampling rate
  
  {//log_lambda: strain sequencing rate parameter
  matrix[N,V] alpha_present;
  for(n in 1:N){
    for(v in 1:V){
      alpha_present[n,v] = alpha[v][2,n];
    }
  }
  for(n in 1:N){
    int i = map_to_location[n];
    int t = map_to_time[n];
    int I = cases[n] > 0 ? 1 : 0;
    real alpha_tot = sum(alpha_present[n,]);
      
    for(v in 1:V){
      real seq_rate;
      real alpha_std = 1/(alpha_tot - alpha_present[n,v] + 1);
      seq_rate = alpha_std*cases[n]+1e-3;
        
      log_lambda[n,v] ~ normal(log(seq_rate), sd_lambda);
        
    }
  }
  }//log_lambda
  
  {//pi: initial presence probability
  
  vector[V] y1_tot;
  for(v in 1:V){
    y1_tot[v] = 0;
    for(i in 1:M){
      y1_tot[v] += y[map_to_n[i,1],v];
    }
  }
  // update to map over 1:M
  for(i in 1:M){ 
    for(v in 1:V) {
      
      int y1 = y[map_to_n[i,1],v];
      real intro = intro_times0[i,v];

      if(y1>0){
        // if observed, high probability of initial presence
        pi_0[i,v] ~ beta(40,1);
      }else if(y1_tot[v]>0){
        // if observed elsewhere, moderate probability of initial presence
        pi_0[i,v] ~ beta(3,4);
      }else if(intro<=0){
        // if already introduced, low/moderate probability of initial presence
        pi_0[i,v] ~ beta(2,10);
      }else{
        // if not observed anywhere, low probability of initial presence
        pi_0[i,v] ~ beta(1,40);
      }
    }
  }
  }//pi
  

  // Spatial weights
  for (n in 1:N_edges) {
    target += log_mix(p_w,
                    normal_lpdf(w[n] | mu_w[1], sd_w[1]),
                    normal_lpdf(w[n] | mu_w[2], sd_w[2]));
  }


  // persistence of a strain ~~> will we want to update this to be country-specific?
  // for(i in 1:M){
  //   log_delta[i] ~ normal(mu_delta, sd_delta);
  // }
  log_delta ~ normal(mu_delta, sd_delta);

  // Gravity parameters
  kappa ~ normal(mu_kappa, sd_kappa); //gravity constant
  zeta_std ~ std_normal();  // distance parameter
    //implies zeta ~ normal(mu_zeta, sd_zeta);  // distance parameter

  for(k in 1:2){
    tau[k] ~ beta(mu_tau[k]*sd_tau, (1-mu_tau[k])*sd_tau); // population parameter
  }
  
  {//gamma_loc: location of introduction
  for(v in 1:N_gamma){
    for(i in 1:N_loc_gamma){
      gamma_loc[v][i] ~ beta(p_gamma_loc[v][i]*l_gamma, (1-p_gamma_loc[v][i])*l_gamma);
    }
  }
  }//gammma_loc
  
  {//gamma_yr: year of introduction
  for(v in 1:N_gamma){
    for(t in 1:N_yr_gamma){
      gamma_yr[v][t] ~ beta(p_gamma_yr[v][t]*l_gamma, (1-p_gamma_yr[v][t])*l_gamma);
      // gamma_yr[v] ~ dirichlet(p_gamma_yr[v]*l_gamma);
    }
  }
  }//gammma_yr

}

generated quantities {
  array[V] matrix<lower = 0, upper = 1>[K,N] beta;    // backward
  matrix<lower = 0, upper = 1>[N,V] rho;    // forward-backward
  matrix<lower =0>[N,V] prev_lambda;    // prevalence accounting for overall probability of presence
  array[N,V] int<lower=1, upper=K> zstar;
  array[N,V] int<lower =0> pred_cases;    // number of cases per strain

  {// start BACKWARD algorithm
  // loop through strains
  for(v in 1:V){

    //initialize
    matrix[K,N] logbeta;    // log backward

    //initialize backward probability at time T
    for(i in 1:M){//start country loop
      int n_T = map_to_n[i,T];
      logbeta[1, n_T] = 0;
      logbeta[2, n_T] = 0;
    }//end country loop

    // start country loop
    for(i in 1:M){
      // start time loop
      for (tforward in 0:(T-2)) {
        
        // mappings
        int tt = T - tforward;
        int n = map_to_n[i,tt]; //n at t for country i (inverted loop)
        int nm = map_to_n[i,tt-1]; //n at t-1 for country i (inverted loop)
        
        // subset to observations in for strain v and country i
        real log_lambda_iv = log_lambda[map_to_n[i,tt],v]; // sequencing rate parameter for poisson approximation of prevalence
        int y_iv = y[map_to_n[i,tt],v];

        // initialize
        real ll0;
        real ll1;
        real a01;
        real a11;
        real a00;
        real a10;

        // transitions from time tt-1 to tt
        a01 = p_col[1][n, v]; // 0 -> 1
        a11 = p_col[2][n, v]; // 1 -> 1

        if (a01 == 0) a01 = epsilon;
        else if (a01 == 1) a01 = 1 - epsilon;

        if (a11 == 0) a11 = epsilon;
        else if (a11 == 1) a11 = 1 - epsilon;

        a00 = 1 - a01; // 0 -> 0
        a10 = 1 - a11; // 1 -> 0
        
        // observation log-probabilities at time tt
        if (y_iv == 0) {
          ll0 = log1m(epsilon);
        } else {
          ll0 = log(epsilon);
        }
        ll1 = poisson_log_lpmf(y_iv | log_lambda[n, v]);
        
        // beta at time tt-1 for state 0
        logbeta[1, nm] = log_sum_exp(logbeta[1, n] + log(a00) + ll0, logbeta[2, n] + log(a01) + ll1);
        // beta at time tt-1 for state 1
        logbeta[2, nm] = log_sum_exp(logbeta[1, n] + log(a10) + ll0,logbeta[2, n] + log(a11) + ll1);

      }// end reverse time loop


      for(t in 1:T){
        // mapping
        int n = map_to_n[i,t];

        // convert beta
        real norm_const = log_sum_exp(logbeta[1, n], logbeta[2, n]);
        beta[v][1, n] = exp(logbeta[1, n] - norm_const);
        beta[v][2, n] = exp(logbeta[2, n] - norm_const);
      }
    }// end country loop

  }// end strain loop

  }// end BACKWARD algorithm


  {//FORWARD-BACKWARD

  // loop through strains
  for(v in 1:V){
    for(t in 1:T){//start time loop
      for(i in 1:M){//start country loop
        // mapping
        int n = map_to_n[i,t];
        real rho0;
        real rho1;
        real denom;

        rho0 = alpha[v][1, n] * beta[v][1, n];
        rho1 = alpha[v][2, n] * beta[v][2, n];
        denom = rho0 + rho1 + 2 * epsilon;
        
        // normalize
        rho[n, v] = (rho1 + epsilon) / denom;

      }//end country loop
    }//end time loop
  }//end strain loop

  }//end forward-backward

  {// prevalence
  matrix[N,V] alphaK;
  for(n in 1:N){

    real rho_0 = prob_v0_calc(n,epsilon,V,rho[n,]); //probability of any strains being present

    if(cases[n]>=1){
        pred_cases[n,] = multinomial_rng(to_vector(lambda_star[n,]), cases[n]);
    }else{
      for(v in 1:V){
        pred_cases[n,v] = 0;
      }

    }

    prev_lambda[n,] = lambda_star[n,]*(1-rho_0);
  }
  }// prevalence

  { // Viterbi algorithm 
  for(i in 1:M){//start country loop
    for(v in 1:V){//start strain loop
      real logp_zstar;
      array[T, 2] int bpointer;         // backpointer to the most likely previous state on the most probable path
      array[T] int zstar_iv;         // z for strain v in i
      array[T, 2] real p;               // max prob for the seq up to t
                                  // with final output from state k for time t

      // t = 1
      {
        int n1 = map_to_n[i, 1];
        int y1 = y[n1, v];
        real ll0;
        real ll1;

        if (y1 == 0) {
          ll0 = log1m(epsilon);
        } else {
          ll0 = log(epsilon);
        }
        ll1 = poisson_log_lpmf(y1 | log_lambda[n1, v]);

        p[1, 1] = log1m(pi_0[i, v]) + ll0;   // state 0
        p[1, 2] = log(pi_0[i, v]) + ll1;     // state 1
      }


      for (t in 2:T) {//start time loop
        int n = map_to_n[i,t]; // n at t for country i
        int y_iv = y[n,v]; // subset observations to country/strain

        // initialize
        real ll0;
        real ll1;
        real a01;
        real a11;
        real a00;
        real a10;

        // transitions 
        a01 = p_col[1][n, v]; // 0 -> 1
        a11 = p_col[2][n, v]; // 1 -> 1

        if (a01 == 0) a01 = epsilon;
        else if (a01 == 1) a01 = 1 - epsilon;

        if (a11 == 0) a11 = epsilon;
        else if (a11 == 1) a11 = 1 - epsilon;

        a00 = 1 - a01; // 0 -> 0
        a10 = 1 - a11; // 1 -> 0
        
        // observation log-probabilities 
        if (y_iv == 0) {
          ll0 = log1m(epsilon);
        } else {
          ll0 = log(epsilon);
        }
        ll1 = poisson_log_lpmf(y_iv | log_lambda[n, v]);

        // current state 0
        {
          real p1 = p[t - 1, 1] + log(a00) + ll0;
          real p2 = p[t - 1, 2] + log(a10) + ll0;
          
          if(p1>p2){
            p[t, 1] = p1;
            bpointer[t,1] = 1;
          }else{
            p[t,1] = p2;
            bpointer[t,1] = 2;
          }
        }
        // current state 1
        {
          real p1 = p[t-1,1] + log(a01) + ll1;
          real p2 = p[t-1,2] + log(a11) + ll1;
          
          if(p1>p2){
            p[t,2] = p1;
            bpointer[t,2] = 1;
          }else{
            p[t,2] = p2;
            bpointer[t,2] = 2;
          }
        }        
      }//end time loop

    logp_zstar = fmax(p[T,1],p[T,2]);

    if (p[T,1] == logp_zstar){
      zstar_iv[T] = 1;
    }else{
      zstar_iv[T] = 2;
    }

    for (t in 1:(T - 1)) {

      zstar_iv[T - t] = bpointer[T - t + 1, zstar_iv[T - t + 1]];
    }

    zstar[map_to_n[i,],v] = zstar_iv;
    }//end strain loop
  }//end country loop

  }//Viterbi

}
