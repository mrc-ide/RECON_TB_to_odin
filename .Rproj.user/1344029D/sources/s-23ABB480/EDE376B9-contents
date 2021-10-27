#Initial states
initial(U) <- N - start_inf
initial(L) <- 0
initial(I) <- start_inf
initial(R) <- 0

#State equations
deriv(U) <- births - U * (lambda + mu)
deriv(L) <- U * lambda * (1-fast) + R * (lambda * (1-fast) * imm) - L * (mu + break_in) 
deriv(I) <- U * lambda * fast + R * (lambda * fast * imm) +  L * break_in + R * relapse - I * (mu + mu_tb + selfcure + Tx)
deriv(R) <- I * (Tx + selfcure) - R * (lambda * imm + relapse + mu)  

# Model Parameters
T_lfx    <- 72                # Life expectancy
T_dur    <- 3                 # Duration of infectious period
beta     <- user(5)                 # Transmission rate per capita
break_in <- 0.1*(1/T_lfx)     # Transition rate from latent into active disease
selfcure <- 0.5*(1/T_dur)     # Rate of spontaneous cure
mu       <- 1/T_lfx           # Background mortality rate
mu_tb    <- 0.5*(1/T_dur)     # TB mortality rate
fast     <- 0.1               # Fraction fast progressing to active TB
imm      <- 0.5               # Infectiousness decline (partial immunity)
relapse  <- 0.005             # Relapse rate
lambda <- beta * I/N          # The force of infection
N <- user(100000)             # Total population
start_inf <- user(1)          # Start infected population
births <- I * mu_tb + N * mu  # Births (for a stable population)
t_intervention <- user(100)   # Intervention time
time_implementation <- user(3)
#This is to simulate a staggered introduction over the time specified in time_implementation
implementation_multiplier <- if(t >= (t_intervention + time_implementation)) 1 else (t - (t_intervention - 1))/time_implementation

#Treatment
T_seek <- user(0)
T_delay <- user(1)

T_cs  <- T_delay * (1 - T_seek)     # Time delay (yrs) between developing symptoms and seeking care
pDx   <- user(0.95)   # Probability of being diagnosed once care sought
pTx   <- user(0.95)   # probability of receiving correct Tx if diagnosed
T_rTx <- user(0.5)    # 6 months treatment duration

Tx <- if(t >= t_intervention) (pDx * pTx * (1/(T_cs + T_rTx))) * implementation_multiplier else 0

#Outputs
output(Incidence)  <- U * (lambda * fast) + R * (lambda * fast * imm) + L * break_in + R * relapse 
output(dIrecent)   <- U * (lambda * fast) + R * (lambda * fast * imm) +  R * relapse                                               
output(dIremote)   <- L * break_in +  R * relapse      