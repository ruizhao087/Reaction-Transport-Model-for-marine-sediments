#################################################################################
# THIS CODE IS DISTRIBUTED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS 
# OF ANY KIND - WITHOUT EVEN THE IMPLIED WARRANTY OF MERCHANTABILITY OR FITNESS
# FOR A PARTICULAR PURPOSE.
#################################################################################

#################################################################################
# PURPOSE: Model the biogeochemical cycling of carbon, nitrogen, oxygen, and Mn 
# species in marine sediment cores from the Arctic Mid-Ocean Ridge.
#
# AUTHORS: Jose M. Mogollon, Rui Zhao
#
#################################################################################

#################################################################################
# Units used: Concentration = mmol; Time = year; Depth = cm.
#################################################################################
# Start model run
# Load R packages
library(ReacTran)
library(marelac)
library(readxl)
library(oce)

# Set work directory 
# setwd("D:/AMOR_sediments_RTM/GS14_GC08_RTM")

#=============================================================================
# Model domain and grid definition
#=============================================================================
L <- 500 # depth of sediment domain [in cm]
N <- 500 # number of grid layers
grid <- setup.grid.1D(x.up = 0, L = L, N = N)  #standard even-spaced grid with N number of nodes

#grid <- setup.grid.1D(x.up = 0, L = L, N = N, dx.1 = 1, p.dx.1 = 1.5)   #Grid setup, see R docs

tmax <- 100000    # number of years for run (only for transient state)
tint <- 100       #time step stored to memory (only for transient state)
tintplot <- 10000 #time step plotted to file (only for transient state)


steadystate <-1      #USE 1 if you want to solve the system under steady state (doesn't always work, but since the sed rates are so slow you would have to run for a very long time in order to run the entire column)
#=============================================================================
# Model parameters:
#=============================================================================
por    <- 0.8 # porosity
porinf <- 0.55 #porosity at infinity
att <- 0.01
svf <- 1.-por # solid volume fraction
S  <- 35 # salinity
TC <- 0 # temperature [deg C]
TK <- TC+273.15   #temperature [K]
P0 <- 1.013 # pressure [bar]
WD <- 2476 #water depth (m)
wdens <- 1.013  #water density (g cm-3)
dwdens <- 2.5    #dry density (g cm-3)
v   <- 0.002   #sedimentation rate (cm yr-1)
            ## Empirical estimation of sedimentation rate (3.3*10^(-0.875-0.000435*WD), cm yr-1, (Middelburg et al, 1996)) are not appliable in our sites.##


Db0 <- 0#2.0     #maximum bioturbation coefficient (usually between 0-10 cm2 yr-1)
x1 <- 3.      #depth of biomixing break (cm)
x2 <-3.       #attenuation of biomixing break (cm)
irr0 <- 0#1.    #maximum bioirrigation coefficient (yr-1)
#=============================================================================
# Pressure Function:
#=============================================================================
Pfunc <- function(x, P0, WD, wdens)
{P <- P0 + (WD*100*980*wdens + x*980*wdens)*1e-6
return(P)}

#=============================================================================
# Porosity Function:
#=============================================================================

Porofunc <- function(x, por, porinf, att, a1, b1) 
{Poro <- a1+b1*(porinf + (por-porinf)*exp(-att*x))
return(Poro)

}

#=============================================================================
# Assign porosity (Poro), solid volume fraction (svf) and pressure (P) grids:
#=============================================================================

P.grid <- setup.prop.1D(func=Pfunc,grid=grid,P0=P0,WD=WD,wdens=wdens)

svf.grid <- setup.prop.1D(func=Porofunc,grid=grid, por=por,porinf=porinf, att=att,a1=1,b1=-1)
Poro.grid <- setup.prop.1D(func=Porofunc,grid=grid,por=por,porinf=porinf,att=att,a1=0,b1=1)


#=============================================================================
# Assign solid and aqueous phase grids: setup.compaction.1D binds Uinf=Vinf, and cannot be used for imposed advective flow (e.g. q=porinf(Uinf-Vinf))
#=============================================================================

Compact.grid <- setup.compaction.1D(v.0=v,por.0=por,por.inf=porinf,por.grid=Poro.grid)

u.grid <- Compact.grid$u
v.grid <- Compact.grid$v


#=============================================================================
# Tortuosity Function:
#=============================================================================

Difffunc <- function(x, Diff, por, porinf, a2, b2, att) {
  poro <- Porofunc(x, por, porinf, att, 0, 1)
  Tortu <- Diff/(a2-b2*log(poro))
  return(Tortu)
}

#=============================================================================
# S function used in biomixing
#=============================================================================
sfunc <- function(x, aval, x1, x2)  {
  return (aval*(exp((x1-x)/x2)/(1.+exp((x1-x)/x2))))
}


#=============================================================================
# Diffusion Coefficients:
#=============================================================================

sinyr <- 60*60*24*365 # number of seconds in one year


Dmol.O2  <- diffcoeff(S = S, t = TC, P = P0)$O2 * sinyr * 10^4  #molecular diffusion coefficient O2
Dmol.Mn  <- diffcoeff(S = S, t = TC, P = P0)$Mn * sinyr * 10^4  #molecular diffusion coefficient Mn
Dmol.NO3 <- diffcoeff(S = S, t = TC, P = P0)$NO3 * sinyr * 10^4  #molecular diffusion coefficient NO3
Dmol.NO2 <- diffcoeff(S = S, t = TC, P = P0)$NO2 * sinyr * 10^4  #molecular diffusion coefficient NO2
Dmol.NH4 <- diffcoeff(S = S, t = TC, P = P0)$NH4 * sinyr * 10^4  #molecular diffusion coefficient NH4
Dmol.HCO3<- diffcoeff(S = S, t = TC, P = P0)$HCO3 * sinyr * 10^4 #molecular diffusion coefficient HCO3
Dmol.SO4 <- diffcoeff(S = S, t = TC, P = P0)$SO4 * sinyr * 10^4  #molecular diffusion coefficient SO4
Dmol.H2S <- diffcoeff(S = S, t = TC, P = P0)$H2S * sinyr * 10^4  #molecular diffusion coefficient H2S
Dmol.PO4 <- diffcoeff(S = S, t = TC, P = P0)$PO4 * sinyr * 10^4  #molecular diffusion coefficient PO4


Db.grid <- setup.prop.1D(func=sfunc, grid=grid, aval=Db0, x1=x1, x2=x2)
alpirr.grid <- setup.prop.1D(func=sfunc, grid=grid, aval=irr0, x1=x1, x2=x2)



# Attachment of parameters to grid
DO2.grid <- setup.prop.1D(func=Difffunc, grid=grid, Diff=Dmol.O2, por=por, porinf=porinf,a2=1, b2=2, att=att)
DMn.grid <- setup.prop.1D(func=Difffunc,grid=grid, Diff=Dmol.Mn, por=por,porinf=porinf,a2=1,b2=2,att=att)
DNO3.grid <- setup.prop.1D(func=Difffunc,grid=grid, Diff=Dmol.NO3, por=por,porinf=porinf,a2=1,b2=2,att=att)
DNO2.grid <- setup.prop.1D(func=Difffunc,grid=grid, Diff=Dmol.NO2, por=por,porinf=porinf,a2=1,b2=2,att=att)
DNH4.grid <- setup.prop.1D(func=Difffunc,grid=grid, Diff=Dmol.NH4, por=por,porinf=porinf,a2=1,b2=2,att=att)
DHCO3.grid <- setup.prop.1D(func=Difffunc,grid=grid, Diff=Dmol.HCO3, por=por,porinf=porinf,a2=1,b2=2,att=att)
DSO4.grid <- setup.prop.1D(func=Difffunc,grid=grid, Diff=Dmol.SO4, por=por,porinf=porinf,a2=1,b2=2,att=att)
DH2S.grid <- setup.prop.1D(func=Difffunc,grid=grid, Diff=Dmol.H2S, por=por,porinf=porinf,a2=1,b2=2,att=att)
DPO4.grid <- setup.prop.1D(func=Difffunc,grid=grid, Diff=Dmol.PO4, por=por,porinf=porinf,a2=1,b2=2,att=att)

#=============================================================================
# Biogeochemical reaction parameters:
#=============================================================================
# 1st order organic matter degradation [1/yr]
  kfox  <- 0.000069  ### Empirical prediction k1=0.4*V^0.6 OR k1=2.2e-5*(FCorg^2.1) (Boudreau book 1997)
  kfox2 <- 0.000005  ### Empirical prediction k2=0.04*V^2 (Boudreau book 1997) 
  kfox3 <- 0.000000  ### most refractory or (0.0=non reactive)
  
# Aerobic respiration inhibition constant [mM]
  kO2 <- 0.005 # h1

# Manganese oxide reduction inhibitation constant [mM]
  kno3 <- 0.01 # h2

  kMnO2 <- 550 # h3 Manganese oxide reduction limitation constant [mM]
  
# Sulfate reduction inhibitation constant [mM]
  kSO4 <- 16
  
# Mnox_rate constant
  k5 <- 110 # k5

# NO3_Mnox_rate constant
  k2 <- 0 # not relevant

# Mn_anamox_rate <- k3*MnO2*NH4*kO2/(O2+kO2)* kno3/(kno3 + NO3)
  k3 <- 0 # not relevant

# nitrification_rate <- k4*NH4*O2
  k4 <- 150 # 

# Anammox_rate constant <- k6*NO3*NH4*kO2/(O2+kO2)
  k6 <- 50 #

# Sulfate reduction decelaration constant
  k8 <- 0.35

#=============================================================================
# Boundary conditions
#=============================================================================
  FCorg   <- 0.0061     # Labile organic matter flux in mol m-2 yr-1 (3. BC)
  FCorg2  <- 0.0032     # refractory organic matter flux in mol m-2 yr-1 (3. BC)
  FCorgNR <- 0.001      # Non-reactive organic matter flux in mol m-2 yr-1 (3. BC)
  FMnO2   <- 0.000002   # flux of manganese oxides mol m-2 yr-1 (3. BC)
  

  O2bw   <- 0.165          # Bottom water oxygen concentration in mM (1. BC)
  Mnbw   <- 0.0001         # Bottom water Mn concentration in mM (1. BC)
  NO3bw  <- 0.025          # Bottom water NO3 concentration in mM (1. BC)
  NO2bw  <- 0.000          # Bottom water NO2 concentration in mM (1. BC)
  NH4bw  <- 0.0001         # Bottom water NH4 concentration in mM (1. BC)
  HCO3bw <- 2.5            # Bottom water HCO3 concentration in mM (1. BC)
  SO4bw  <- 29.0           # Bottom water SO4 concentration in mM (1. BC)
  H2Sbw  <- 0.0001         # Bottom water H2S concentration in mM (1. BC)
  PO4bw  <- 0.0026         # Bottom water PO4 concentration in mM (1. BC)

# Flux conversion to umol cm-2 to ultimately keep all concentrations consistent in mmol L-1
  FCorg  <- FCorg * 100
  FCorg2 <- FCorg2 * 100
  FMnO2  <- FMnO2 * 100

  FCorgNR <- FCorgNR*100.



#=============================================================================
# Model formulation
#=============================================================================

AMORmodel <- function(t, state, pars) {
  # Initialisation of state variables
  
  Corg     <-  state[(0*N+1):(1*N)]
  O2       <-  state[(1*N+1):(2*N)]
  Mn       <-  state[(2*N+1):(3*N)]
  NO3      <-  state[(3*N+1):(4*N)]
  MnO2     <-  state[(4*N+1):(5*N)]
  NH4      <-  state[(5*N+1):(6*N)]
  Corg2    <-  state[(6*N+1):(7*N)]
  HCO3     <-  state[(7*N+1):(8*N)]
  SO4      <-  state[(8*N+1):(9*N)]
  PO4      <-  state[(9*N+1):(10*N)]
  H2S      <-  state[(10*N+1):(11*N)]
  CorgNR   <-  state[(11*N+1):(12*N)]
  NO2      <-  state[(12*N+1):(13*N)]
  
  #### Updating time varying variables #############################################################
  # Nothing to update, SS (steady state) model
  
  
  #### Define transport component #############################################################
  tranCorg <- tran.1D(C = Corg, flux.up = FCorg, D = Db.grid, v = v.grid,
                      VF = svf.grid, dx = grid)
  tranO2   <- tran.1D(C = O2, C.up = O2bw, D = DO2.grid, v = u.grid,
                     VF = Poro.grid, dx = grid)
  tranMn   <- tran.1D(C = Mn, C.up = Mnbw, D = DMn.grid, v = u.grid,
                    VF = Poro.grid, dx = grid)
  tranNO3  <- tran.1D(C = NO3, C.up = NO3bw, D = DNO3.grid, v = u.grid,
                     VF = Poro.grid, dx = grid)
  tranMnO2 <- tran.1D(C = MnO2, flux.up = FMnO2, D = Db.grid, v = v.grid,
                      VF = svf.grid, dx = grid)
  tranNH4  <- tran.1D(C = NH4, C.up = NH4bw, D = DNH4.grid, v = u.grid,
                     VF = Poro.grid, dx = grid)
  tranCorg2<- tran.1D(C = Corg2, flux.up = FCorg2, D = Db.grid, v = v.grid,
                       VF = svf.grid, dx = grid)
  tranHCO3 <- tran.1D(C = HCO3, C.up = HCO3bw, D = DHCO3.grid, v = u.grid,
                      VF = Poro.grid, dx = grid)
  tranSO4  <- tran.1D(C = SO4, C.up = SO4bw, D = DSO4.grid, v = u.grid,
                     VF = Poro.grid, dx = grid)
  tranPO4  <- tran.1D(C = PO4, C.up = PO4bw, D = DPO4.grid, v = u.grid,
                     VF = Poro.grid, dx = grid)
  tranH2S  <- tran.1D(C = H2S, C.up = H2Sbw, D = DH2S.grid, v = u.grid,
                     VF = Poro.grid, dx = grid)
  tranCorgNR <- tran.1D(C = CorgNR, flux.up = FCorgNR, D = Db.grid, v = v.grid,
                        VF = svf.grid, dx = grid)
  
  #### Define reaction rates #############################################################
  C_org_rate  <- kfox*Corg
  C_org_rate2 <- kfox2*Corg2
  C_org_rate3 <- kfox3*CorgNR
  
  O2_Coxrate1   <- (C_org_rate) *O2 / (O2+kO2)
  NO3_Coxrate1  <- (C_org_rate) *kO2 / (O2+kO2) * NO3/(kno3 + NO3)
  MnO2_Coxrate1 <- (C_org_rate) *kO2 / (O2+kO2) * kno3/(kno3 + NO3) *MnO2/(kMnO2 + MnO2)
  SO4_Coxrate1  <- k8*(C_org_rate) *kO2 / (O2+kO2) * kno3/(kno3 + NO3) * SO4/(kSO4+SO4)
  
  O2_Coxrate2   <- (C_org_rate2) *O2 / (O2+kO2)
  NO3_Coxrate2  <- (C_org_rate2) *kO2 / (O2+kO2) * NO3/(kno3 + NO3)
  MnO2_Coxrate2 <- (C_org_rate2) *kO2 / (O2+kO2) * kno3/(kno3 + NO3) *MnO2/(kMnO2 + MnO2)
  SO4_Coxrate2  <- k8*(C_org_rate2) *kO2 / (O2+kO2) * kno3/(kno3 + NO3) * SO4/(kSO4+SO4) 
  
  O2_Coxrate3   <- (C_org_rate3) *O2 / (O2+kO2)
  NO3_Coxrate3  <- (C_org_rate3) *kO2 / (O2+kO2) * NO3/(kno3 + NO3)
  MnO2_Coxrate3 <- (C_org_rate3) *kO2 / (O2+kO2) * kno3/(kno3 + NO3) *MnO2/(kMnO2 + MnO2)
  SO4_Coxrate3  <- k8*(C_org_rate3) *kO2 / (O2+kO2) * kno3/(kno3 + NO3) * SO4/(kSO4+SO4) 
  
  O2_Coxrate   <- O2_Coxrate1 + O2_Coxrate2+ O2_Coxrate3
  NO3_Coxrate  <- NO3_Coxrate1 + NO3_Coxrate2+ NO3_Coxrate3
  MnO2_Coxrate <- MnO2_Coxrate1 + MnO2_Coxrate2+ MnO2_Coxrate3
  SO4_Coxrate  <- SO4_Coxrate1+ SO4_Coxrate2+ SO4_Coxrate3
  
  Mnox_rate      <- k5*Mn*O2
  NO3_Mnox_rate  <- k2*Mn*NO3*kO2/(O2+kO2)
  Mn_anamox_rate <- k3*MnO2*NH4*kO2/(O2+kO2)* kno3/(kno3 + NO3)
  nitrification_rate <- k4*NH4*O2
  Anammox_rate   <- k6*NO3*NH4*kO2/(O2+kO2)
  
  
  #### Define derivatives #############################################################
  q <- svf.grid$mid / Poro.grid$mid #1/svf to 1/por
  r <- Poro.grid$mid / svf.grid$mid #1/por to 1/svf
  
  dCorgdt <- tranCorg$dC - O2_Coxrate1 - NO3_Coxrate1 - MnO2_Coxrate1 -SO4_Coxrate1
  
  dCorg2dt <- tranCorg2$dC - O2_Coxrate2 - NO3_Coxrate2 - MnO2_Coxrate2 -SO4_Coxrate2
  
  dCorgNRdt <- tranCorgNR$dC - O2_Coxrate3 - NO3_Coxrate3 - MnO2_Coxrate3 -SO4_Coxrate3

  dO2dt  <- tranO2$dC - q*O2_Coxrate - Mnox_rate - 2*nitrification_rate +(O2bw-O2)*alpirr.grid$mid
  
  dMndt  <- tranMn$dC + q*2.*MnO2_Coxrate - 2*Mnox_rate - NO3_Mnox_rate + Mn_anamox_rate + (Mnbw-Mn)*alpirr.grid$mid
  
  dNO3dt  <- tranNO3$dC - q*4./5.*NO3_Coxrate + nitrification_rate - 2./5.*NO3_Mnox_rate -3./5.*Anammox_rate + (NO3bw-NO3)*alpirr.grid$mid
  
  dMnO2 <- tranMnO2$dC - 2.0*MnO2_Coxrate + 2*r*Mnox_rate + r*NO3_Mnox_rate - r*Mn_anamox_rate
  
  dNH4dt <- tranNH4$dC + q*1./5.*(O2_Coxrate +NO3_Coxrate+ MnO2_Coxrate+ SO4_Coxrate) - nitrification_rate  -2./3.*Mn_anamox_rate -Anammox_rate +(NH4bw-NH4)*alpirr.grid$mid
  
  dHCO3dt <-tranHCO3$dC + q*(O2_Coxrate +NO3_Coxrate+MnO2_Coxrate +SO4_Coxrate) +(HCO3bw-HCO3)*alpirr.grid$mid
  
  dSO4dt <-tranSO4$dC - q*0.5*SO4_Coxrate +(SO4bw-SO4)*alpirr.grid$mid 
  
  dPO4dt <-tranPO4$dC + q*1/106.*(O2_Coxrate +MnO2_Coxrate+ +NO3_Coxrate+ SO4_Coxrate) +(PO4bw-PO4)*alpirr.grid$mid 
  
  dH2Sdt <-tranH2S$dC + q*0.5*SO4_Coxrate +(H2Sbw-H2S)*alpirr.grid$mid 
  
 
  #### Return variables ################################################################
  return(list(c(dCorgdt=dCorgdt,dO2dt=dO2dt, dMndt=dMndt, dNO3dt=dNO3dt, dMnO2=dMnO2, 
                dNH4dt=dNH4dt,dCorg2dt=dCorg2dt,dHCO3dt=dHCO3dt,dSO4dt=dSO4dt, dPO4dt=dPO4dt, 
                dH2Sdt=dH2Sdt, dCorgNRdt=dCorgNRdt), O2_Coxrate=O2_Coxrate, NO3_Coxrate=NO3_Coxrate,
              MnO2_Coxrate=MnO2_Coxrate,SO4_Coxrate=SO4_Coxrate,Mnox_rate=Mnox_rate,
              NO3_Mnox_rate=NO3_Mnox_rate,Mn_anamox_rate=Mn_anamox_rate,
              nitrification_rate=nitrification_rate, Anammox_rate=Anammox_rate))
  
} ## end of AMORmodel function



########################################
# Initial conditions of state variables:
########################################
Corg.in <- rep(1e-6, length.out = N)
O2.in   <- rep(1e-7, length.out = N)
Mn.in   <- rep(1e-6, length.out = N)
NO3.in  <- rep(0.03, length.out = N)
MnO2.in <- rep(1e-10, length.out = N)
NH4.in <- rep(0.0001, length.out = N)
Corg2.in <- rep(1e-8, length.out = N)
HCO3.in <- rep(0.1, length.out = N)
SO4.in <- rep(27., length.out = N)
PO4.in <- rep(0.0026, length.out = N)
H2S.in <- rep(1e-8, length.out = N)
CorgNR.in <-  rep(1e-6, length.out = N)

state <- c(Corg.in, O2.in, Mn.in, NO3.in, MnO2.in, NH4.in,Corg2.in,HCO3.in,SO4.in, PO4.in, H2S.in, CorgNR.in)
names <- c("Corg", "O2", "Mn", "NO3", "MnO2", "NH4","Corg2","HCO3","SO4", "PO4", "H2S","CorgNR")


if(steadystate==0) {
  ########################################
  # Transient model run
  ########################################
  times <- seq(0, tmax , by = tint)
  print(system.time(
    std <- ode.1D(y = state, times = times, func = AMORmodel, parms = NULL,
                  names = names, method = "lsoda", verbose = TRUE, nspec = length(names), atol= 1e-20)
  ))
  
  outLoop <- seq(0,tmax, by = tintplot)
  for(i in 2:length(outLoop)){
    row    <- floor(outLoop[i]/tint)
    
    
    
    Corg      <- std[row, (0*N+2):(1*N+1)]
    O2        <- std[row, (1*N+2):(2*N+1)]
    Mn        <- std[row, (2*N+2):(3*N+1)]
    NO3       <- std[row, (3*N+2):(4*N+1)]
    MnO2      <- std[row, (4*N+2):(5*N+1)]
    NH4       <- std[row, (5*N+2):(6*N+1)]
    Corg2     <- std[row, (6*N+2):(7*N+1)]
    
    O2_Coxrate     <- std[row, (7*N+2):(8*N+1)]
    NO3_Coxrate    <- std[row, (8*N+2):(9*N+1)]
    MnO2_Coxrate   <- std[row, (9*N+2):(10*N+1)]
    Mnox_rate      <- std[row, (10*N+2):(11*N+1)]
    NO3_Mnox_rate  <- std[row, (11*N+2):(12*N+1)]
    Mn_anamox_rate <- std[row, (12*N+2):(13*N+1)]
    nitrification_rate <- std[row, (13*N+2):(14*N+1)]
    Anammox_rate   <- std[row, (14*N+2):(15*N+1)]
    SO4_Coxrate    <- std[row, (15*N+2):(16*N+1)]
    
    
    depth <- grid$x.mid
    pdf(paste(c("test_new","_",row*tint,".pdf"), collapse = ""), width = 15)
    par(mfrow=c(3,4), mar=c(4, 2, 1, 1))
    
    # calculate depth-integrated rates, convert from mol e- m-2 s-1 to mol C m-2 s-1
    O2Rint <- integrateTrapezoid(depth,Poro.grid$mid*(O2_Coxrate))
    MNOXint <- integrateTrapezoid(depth,Poro.grid$mid*(MnO2_Coxrate))
    NO3REDint <- integrateTrapezoid(depth,Poro.grid$mid*(NO3_Coxrate))
    SO4REDint <- integrateTrapezoid(depth,Poro.grid$mid*(SO4_Coxrate))
    
    # plot results
    
    plot((Corg+Corg2)*1e-6/2.6*12.01*100, depth, ylim=c(max(depth), min(depth)), type = "l", col="cyan",
         xlab="Corg (wt%)", ylab="depth [cm]", xlim = c(0,1.0))
    
    lines((MnO2)*1e-6/2.6*54.938*100./10., depth, ylim=c(max(depth), min(depth)), type = "l", col="red",
          xlab="Mn (wt%)", ylab="depth [cm]", xlim = c(0,0.1))
    
    plot(O2, depth, ylim=c(max(depth), min(depth)), type = "l", col="magenta",
         xlab="O2 (mM)", ylab="depth [cm]", xlim = c(0,0.200))
    
    plot(Mn, depth, ylim=c(max(depth), min(depth)), type="l", col="blue",
         xlab="Mn (mM)", ylab="depth [cm]", xlim = c(0,3.0))
    
    plot(NO3, depth, ylim=c(max(depth), min(depth)), type="l", col="green",
         xlab="NO3 (mM)", ylab="depth [cm]", xlim = c(0,0.1))
    
    plot(NH4, depth, ylim=c(max(depth), min(depth)), type="l", col="green",
         xlab="NH4 (mM)", ylab="depth [cm]", xlim = c(0,1.0))
    
    plot(O2_Coxrate, depth, ylim=c(max(depth), min(depth)), type="l", xlab="O2R (mM C yr-1)", ylab="depth [cm]")
    plot(MnO2_Coxrate, depth, ylim=c(max(depth), min(depth)), type="l", xlab="Mnox (mM C yr-1)", ylab="depth [cm]")
    plot(NO3_Coxrate, depth, ylim=c(max(depth), min(depth)), type="l", xlab="NO3red (mM C yr-1)", ylab="depth [cm]")
    plot(Mnox_rate, depth, ylim=c(max(depth), min(depth)), type="l", xlab="Mn2+ O2-oxidation (mM O2 yr-1)", ylab="depth [cm]")
    plot(NO3_Mnox_rate, depth, ylim=c(max(depth), min(depth)), type="l", xlab="Mn2+ NO3-oxidation (mM Mn yr-1)", ylab="depth [cm]")
    plot(Mn_anamox_rate, depth, ylim=c(max(depth), min(depth)), type="l", xlab="Mn-anammox (mM Mn yr-1)", ylab="depth [cm]")
    plot(nitrification_rate, depth, ylim=c(max(depth), min(depth)), type="l", xlab="Denitrif. (mM N yr-1)", ylab="depth [cm]")
    plot(Anammox_rate, depth, ylim=c(max(depth), min(depth)), type="l", xlab="Anammox. (mM N yr-1)", ylab="depth [cm]")
    
    #plot(Poro.grid$mid, depth, ylim=c(max(depth), min(depth)), type = "l", col="grey",
    #xlab="Porosity", ylab="depth [cm]")
    
    dev.off()
    
    
  }
} else {
  ########################################
  # Steady State
  ########################################
  std <- steady.1D(y = state, func = AMORmodel, parms = NULL,
                   names = names, nspec = length(names), pos=TRUE)
  
  # extract results
  Corg  <- std$y[,1]
  O2    <- std$y[,2]
  Mn    <- std$y[,3]
  NO3   <- std$y[,4]
  MnO2  <- std$y[,5]
  NH4   <- std$y[,6]
  Corg2 <- std$y[,7]
  HCO3  <- std$y[,8]
  SO4   <- std$y[,9]
  PO4   <- std$y[,10]
  H2S   <- std$y[,11]
  CorgNR<- std$y[,12]
  
  O2_Coxrate   <- std$O2_Coxrate
  NO3_Coxrate  <- std$NO3_Coxrate
  MnO2_Coxrate <- std$MnO2_Coxrate
  SO4_Coxrate  <- std$SO4_Coxrate
  
  Mnox_rate    <- std$Mnox_rate
  NO3_Mnox_rate <- std$NO3_Mnox_rate
  Mn_anamox_rate <- std$Mn_anamox_rate
  nitrification_rate <- std$nitrification_rate
  Anammox_rate   <- std$Anammox_rate
  SO4_Coxrate    <- std$SO4_Coxrate
  
  #### Prepare the model output and save as a csv file!  
#  depth <- grid$x.mid
#  Output_data <- std$y
#  Output_data <- cbind(Output_data,O2_Coxrate,NO3_Coxrate,MnO2_Coxrate,Anammox_rate,SO4_Coxrate,nitrification_rate,depth)

  ## Export model for check or plotting ##   
#  write.csv(Output_data, "D:/North Pond Maintanence energy/Nitrifiers_in_oxic_sediments/GC08_RTM/GC08_model_output.csv")
  
  ## Load the measured profiles #### Remember to change it to the exact path on your machine.
  GC08data <- read_excel("D:/AMOR_sediments_RTM/AMOR_profiles.xlsx",
                       sheet="GC08") 
  
  depth <- grid$x.mid
  pdf(paste(c("test_new_ss_GC08.pdf"), collapse = ""), width = 9)
  par(mfrow=c(3,6), mar=c(4, 2, 1, 1))
  
  # calculate depth-integrated rates, convert from mol e- m-2 s-1 to mol C m-2 s-1
  O2Rint <- integrateTrapezoid(depth,Poro.grid$mid*(O2_Coxrate))
  MNOXint <- integrateTrapezoid(depth,Poro.grid$mid*(MnO2_Coxrate))
  NO3REDint <- integrateTrapezoid(depth,Poro.grid$mid*(NO3_Coxrate))
  
  depthy <- 400.0 #max(depth)
  
  depth_OrgC <- GC08data$Depth_CN
  OrgC_vals  <- GC08data$Organic_carbon
  plot((Corg+Corg2+CorgNR)*1e-6/2.6*12.*100, depth, ylim=c(depthy, min(depth)), type = "l",lwd=2, col="cyan",
       xlab="Corg (wt%-cyan), MnO2 (wt%/10-red)", ylab="depth [cm]", xlim = c(0, 2.0))
  points(OrgC_vals, depth_OrgC, pch=16)
  
  lines((MnO2)*1e-6/2.6*54.938*100, depth, ylim=c(depthy, min(depth)), type = "l", col="red",
        xlab="Mn (wt%)", ylab="depth [cm]", xlim = c(0,0.1))

  
  depth_O2 <- GC08data$Depth_O2
  O2_vals  <- GC08data$Oxygen
  plot(O2, depth, ylim=c(depthy, min(depth)), type = "l",lwd=1, col="magenta",
       xlab="O2 (mM)", ylab="depth [cm]", xlim = c(0,0.200))
  points(O2_vals, depth_O2, pch=16, col="magenta")
  
  depth_Mn <- GC08data$Depth
  Mn_vals  <- GC08data$Mn
  plot(Mn, depth, ylim=c(depthy, min(depth)), type = "l", lwd=2,col="blue",
       xlab="Mn (mM)", ylab="depth [cm]", xlim = c(0,0.2))
  points(Mn_vals, depth_Mn, pch=16,col="blue")
  
  depth_NO3 <- GC08data$Depth
  NO3_vals  <- GC08data$Nitrate
  plot(NO3, depth, ylim=c(depthy, min(depth)), type = "l",lwd=2, col="green",
       xlab="NO3 (mM)", ylab="depth [cm]", xlim = c(0,0.04))
  points(NO3_vals, depth_NO3, pch=16,col="green")
  
  depth_NH4 <- GC08data$Depth
  NH4_vals <- GC08data$NH4
  plot(NH4, depth, ylim=c(depthy, min(depth)), type = "l", lwd=2,col="orange",
       xlab="NH4 (mM)", ylab="depth [cm]", xlim = c(0,0.06))
  points(NH4_vals, depth_NH4, pch=16,col="orange")
  
  depth_HCO3 <- GC08data$Depth
  HCO3_vals  <- GC08data$DIC
  plot(HCO3, depth, ylim=c(depthy, min(depth)), type = "l", lwd=2,col="firebrick",
       xlab="DIC (mM)", ylab="depth [cm]", xlim = c(2, 4))
  points(HCO3_vals, depth_HCO3, pch=16,col="firebrick")
  
  depth_alk <- GC08data$Depth
  alk_vals  <- GC08data$Alkalinity
  plot(alk_vals, depth_alk, ylim=c(depthy, min(depth)), col="chocolate",
       xlab="Alkanility (mM)", ylab="depth [cm]", xlim=c(2, 5))
  
  depth_pH <- GC08data$Depth
  pH_vals  <- GC08data$pH
  plot(pH_vals, depth_pH, ylim=c(depthy, min(depth)), col="deeppink",
       xlab="pH", ylab="depth [cm]", xlim=c(7, 8))
  
  depth_Fe <- GC08data$Depth
  Fe_vals  <- GC08data$`Iron(II)`
  plot(Fe_vals, depth_Fe, ylim=c(depthy, min(depth)), col="black",
       xlab="Fe(II) (mM)", ylab="depth [cm]", xlim=c(0, 0.2))
  
  depth_SO4 <- GC08data$Depth
  SO4_vals  <- GC08data$Sulfate
  plot(SO4, depth, ylim=c(depthy, min(depth)), type = "l",lwd=2, col="purple",
       xlab="SO4 (mM)", ylab="depth [cm]", xlim = c(25,30.0))
  points(SO4_vals, depth_SO4, pch=16, col="purple")
  
  
  plot(H2S, depth, ylim=c(depthy, min(depth)), type = "l", lwd=2,col="purple",
       xlab="H2S (mM)", ylab="depth [cm]", xlim = c(0,0.2))
  
  depth_PO4 <- GC08data$Depth
  PO4_vals  <- GC08data$PO4
  plot(PO4, depth, ylim=c(depthy, min(depth)), type = "l", lwd=2,col="chocolate",
       xlab="PO4 (mM)", ylab="depth [cm]", xlim = c(0,0.02))
  points(PO4_vals, depth_PO4)
  
  depth_Si <- GC08data$Depth
  Si_vals  <- GC08data$Si
  plot(Si_vals, depth_Si, ylim=c(depthy, min(depth)), col="black",
       xlab="Si (mM)", ylab="depth [cm]", xlim = c(0, 0.3))
  
  depth_OrgN <- GC08data$Depth_CN
  OrgN_vals  <- GC08data$Organic_Nitrogen
  plot(OrgN_vals, depth_OrgN, xlim=c(0,0.2), ylim=c(depthy, min(depth)), xlab="Org N (wt%)")
  
  CN_ratio  <- GC08data$C_N_ratio
  plot(CN_ratio, depth_OrgN, xlim=c(0,10), ylim=c(depthy, min(depth)), xlab="C/N")
  
  depth_poro <- GC08data$Depth_CN
  poro_vals  <- GC08data$Porosity
  plot(Poro.grid$mid, depth, ylim=c(depthy, min(depth)), type = "l", lwd=2, col="black",
       xlab="Porosity", ylab="depth [cm]", xlim = c(0.4,1.0))
  points(poro_vals, depth_poro)
  
  plot(CorgNR, depth, ylim=c(depthy, min(depth)), type = "l", col="purple",
       xlab="CorgNR (mM)", ylab="depth [cm]", xlim = c(0, 0.01))
  
  ###reaction rates plotting###
  
  plot(O2_Coxrate, depth, ylim=c(depthy, min(depth)), type="l", xlab="O2R (µmol cm-3 C yr-1)", ylab="depth [cm]", col="magenta", lwd=2)
  plot(MnO2_Coxrate, depth, ylim=c(depthy, min(depth)), type="l", xlab="Mnox reduction (µmol cm-3 C yr-1)", ylab="depth [cm]")
  plot(NO3_Coxrate, depth, ylim=c(depthy, min(depth)), type="l", xlab="NO3red (µmol cm-3 C yr-1)", ylab="depth [cm]", col="green", lwd=2)
  plot(Mnox_rate, depth, ylim=c(depthy, min(depth)), type="l", xlab="Mn2+ O2-oxidation (µmol cm-3 O2 yr-1)", ylab="depth [cm]")
  plot(NO3_Mnox_rate, depth, ylim=c(depthy, min(depth)), type="l", xlab="Mn2+ NO3-oxidation (µmol cm-3 Mn yr-1)", ylab="depth [cm]")
  plot(Mn_anamox_rate, depth, ylim=c(depthy, min(depth)), type="l", xlab="Mn-anammox (µmol cm-3 Mn yr-1)", ylab="depth [cm]")
  plot(nitrification_rate, depth, ylim=c(depthy, min(depth)), type="l", xlab="Nitrification (µmol cm-3 N yr-1)", ylab="depth [cm]", col="orange", lwd=2)
  plot(Anammox_rate, depth, ylim=c(depthy, min(depth)), type="l", xlab="Anammox. (µmol cm-3 N yr-1)", ylab="depth [cm]",col="red", lwd=2)
  plot(SO4_Coxrate, depth, ylim=c(depthy, min(depth)), type="l", xlab="Sulfate reduction (µmol cm-3 C yr-1)", ylab="depth [cm]", col="purple", lwd=2) 
  plot(log(NO3_Coxrate/Anammox_rate), depth, ylim=c(200, min(depth)), type="l", xlab="Log (Denitrification/Anammox)", ylab="depth [cm]", col="purple", lwd=2) 
  
  dev.off()
}
