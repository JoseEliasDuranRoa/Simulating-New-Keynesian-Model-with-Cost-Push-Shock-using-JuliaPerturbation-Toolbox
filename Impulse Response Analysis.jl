## José Elías Durán Roa ##
## PART 1 - IMPULSE RESPONSE ANALYSIS ESTIMATION ## 

###################################################################################
## STEP N°1 = First we modify the general equations and description of the model ##
###################################################################################

using MKL
include("solution_functions_7_0.jl");

## Model
    # Adjustment
        flag_order      = 1
        flag_deviation  = true
        flag_SSsolver   = false

    # Parameters
        @vars THETA BETTA GAMMA DELTA ALPHA XI RHO PHI_R PHI_PPI PHI_Y PPS YSS RRS SIGMAA SIGMAETA SIGMAMP RHOETA sseta
        parameters  =   [THETA; BETTA; GAMMA; DELTA; ALPHA; XI; RHO; PHI_R; PHI_PPI; PHI_Y; PPS; YSS; RRS; SIGMAA; SIGMAETA; SIGMAMP; RHOETA ; sseta]
        estimate    =   [RHO; SIGMAA; SIGMAETA; SIGMAMP; RHOETA]
        position = [7; 14; 15; 16; 17]
        priors = (prior_rho = (iv = .9 , lb = 0.0 , ub = 1.0, d = truncated(Normal(0.75, 0.25),0.0,1.0)),
                  prior_sigma = (iv = .05, lb = 0.0, ub = 1e6, d = InverseGamma(5.0, 0.25)),prior_sigmaeta = (iv = .05, lb = 0.0, ub = 1e6, d = InverseGamma(5.0, 0.25)), 
                  prior_sigmamp = (iv = .05, lb = 0.0, ub = 1e6, d = InverseGamma(5.0, 0.25)), 
                  prior_rhoeta = (iv = .9 , lb = 0.0 , ub = 1.0, d = truncated(Normal(0.75, 0.25),0.0,1.0)))
 

    # Variables
        @vars c cp h hp l lp w wp mc mcp s sp ph php ppi ppip rt rtp x1 x1p x2 x2p br brp brb brbp z zp k kp ii iip u up etas etasp
        x = [s; z; brb; k; etas]
        y = [c; h; ppi; rt; br; w; l; mc; ph; x1; x2; ii; u]
        xp = [sp; zp; brbp; kp; etasp]
        yp = [cp; hp; ppip; rtp; brp; wp; lp; mcp; php; x1p; x2p; iip; up]
        variables = [x; y; xp; yp]



    # Shock
        @vars epsilon
        e = [epsilon]
        eta = Array([0.0; SIGMAA ; SIGMAMP ; 0.0 ; SIGMAETA])

    # Equilibrium conditions
        f1 = -kp + (1-DELTA)*k + ii
        f2  = l - (c*(1-h)^GAMMA)^(-XI) * (1-h)^GAMMA
        f3  = BETTA*lp/ppip - l*rt
        f4  = 1/br - rt
        f5  = -sp + (1-THETA)*(ph^(-etas)) + THETA*(ppi^etas)*s
        f6  = z*h^(1-ALPHA)*k^(ALPHA) - sp*(c+ii)
        f7  = -1 + (1-THETA)*(ph^(1-etas)) + THETA*(ppi^(etas-1))
        f8  = -x1 + (ph^(1-etas))*c*((etas-1)/etas) + THETA*rt*((ph/php)^(1-etas))*(ppip^(etas))*x1p
        f9  = -x2 + (ph^(-etas))*c*mc + THETA*rt*((ph/php)^(-etas))*(ppip^(etas+1))*x2p
        f10 = x2-x1
        f11 = w - mc*z*(1-ALPHA)*h^(-ALPHA)*k^(ALPHA)
        f12 = w - ((c*(1-h)^GAMMA)^(-XI)*GAMMA*(1-h)^(GAMMA-1)*c)/((c*(1-h)^GAMMA)^(-XI)*(1-h)^GAMMA)
        f13a = u - mc*z*ALPHA*h^(1-ALPHA)*k^(ALPHA-1)
        f13b = l - BETTA*(up+1-DELTA)*lp
        f14a = log(br/RRS) - PHI_R*log(brb/RRS) - PHI_PPI*log(ppi/PPS) - PHI_Y*log((sp*(c+ii))/YSS)
        f14b = brbp - br
        f15 = log(zp) - RHO*log(z)
        f16 = log(etasp) - (1-RHOETA)*log(sseta)-RHOETA*log(etas)

        f = [f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11;f12;f13a;f13b;f14a;f14b;f15;f16]

        

    # Steady state
        # Values
       
        Z           =       1.00
        sseta       =       3.00
        R           =       1/RRS
        RR          =       RRS
        PH          =       ((1-THETA*PPS^(sseta-1))/(1-THETA))^(1/(1-sseta))
        PPI         =       PPS
        S           =       (1-THETA)*PH^(-sseta)/(1-THETA*PPS^(sseta))
        MC          =       ((1-THETA*RRS^(-1)*PPS^(sseta+1))*((sseta-1)/sseta)*PH)/(1-THETA*RRS^(-1)*PPS^(sseta))
        U           =       1/BETTA - 1 + DELTA
        H           =       ((1-ALPHA)*MC*Z*((ALPHA*Z*MC)/U)^(ALPHA/(1-ALPHA)))/((1-ALPHA)*MC*Z*((ALPHA*Z*MC)/U)^(ALPHA/(1-ALPHA)) + (GAMMA*Z*((ALPHA*Z*MC)/U)^(ALPHA/(1-ALPHA)) - S*DELTA*((ALPHA*Z*MC)/U)^(1/(1-ALPHA)))/S);
        K           =       ((ALPHA*Z*MC)/U)^(1/(1-ALPHA))*H
        W           =       (1-ALPHA)*MC*Z*(K/H)^ALPHA
        Y           =       Z*H^(1-ALPHA)*K^ALPHA
        YSS         =       Y
        C           =       (Y - S*DELTA*K)/S
        II          =       DELTA*K
        X1          =       (PH^(1-sseta)*Y*S^(-1)*((sseta-1)/sseta))/(1-THETA*RRS^(-1)*PPS^sseta)
        X2          =       X1
        LAM         =       C^(-XI)*(1-H)^(GAMMA*(1-XI))
        

        # Vector
        SS = [
                log(S);
                log(Z);
                log(sseta);
                log(RR);
                log(K);
                log(C);
                log(H);
                log(PPI);
                log(R);
                log(RR);
                log(W);
                log(LAM);
                log(MC);
                log(PH);
                log(X1);
                log(X2);
                log(II);
                log(U);
                log(S);
                log(Z);
                log(sseta);
                log(RR);
                log(K);
                log(C);
                log(H);
                log(PPI);
                log(R);
                log(RR);
                log(W);
                log(LAM);
                log(MC);
                log(PH);
                log(X1);
                log(X2);
                log(II);
                log(U);
            ]

        PAR_SS = [THETA; BETTA; GAMMA; DELTA; ALPHA; XI; RHO; PHI_R; PHI_PPI; PHI_Y; PPS; YSS; RRS; SIGMAA; SIGMAETA; SIGMAMP; RHOETA; sseta]


        # Procesing the model (no adjustment needed) 
        model = (parameters = parameters, estimate = estimate, estimation = position,
                    npar = length(parameters), ns = length(estimate), priors = priors,
                    x = x, y = y, xp = xp, yp = yp, variables = variables,
                    nx = length(x), ny = length(y), nvar = length(x) + length(y),
                    e = e, eta = eta,
                    ne = length(e),
                    f = f,
                    nf = length(f),
                    SS = SS, PAR_SS = PAR_SS,
                    flag_order = flag_order, flag_deviation = flag_deviation, flag_SSsolver = flag_SSsolver)

        process_model(model)



## Solution
    # Parametrization
        THETA      =       0.80
        BETTA      =       1.04^(-1/4)
        GAMMA      =       3.6133
        DELTA      =       0.01
        ALPHA      =       0.30
        XI         =       2.00
        RHO        =       0.8556
        Z          =       1.00
        sseta      =       3.00
        PHI_PPI    =       3.00
        PHI_Y      =       0.01
        PHI_R      =       0.80
        PPS        =       1.042^(1/4)
        RRS        =       PPS/BETTA
        SIGMAA     =       0.05
        SIGMAETA   =       0.05
        SIGMAMP    =       0.05
        RHOETA     =       0.825

        
        Z          =       1.00
        sseta      =       3.00
        R          =       1.00/RRS
        RR         =       RRS
        PH         =       ((1.00-THETA*PPS^(sseta-1.00))/(1.00-THETA))^(1.00/(1.00-sseta))
        PPI        =       PPS
        S          =       (1.00-THETA)*PH^(-sseta)/(1.00-THETA*PPS^(sseta))
        MC         =       ((1.00-THETA*RRS^(-1.00)*PPS^(sseta+1.00))*((sseta-1.00)/sseta)*PH)/(1.00-THETA*RRS^(-1.00)*PPS^(sseta))
        U          =       1.00/BETTA - 1.00 + DELTA
        H          =       ((1.00-ALPHA)*MC*Z*((ALPHA*Z*MC)/U)^(ALPHA/(1.00-ALPHA)))/((1.00-ALPHA)*MC*Z*((ALPHA*Z*MC)/U)^(ALPHA/(1.00-ALPHA)) + (GAMMA*Z*((ALPHA*Z*MC)/U)^(ALPHA/(1.00-ALPHA)) - S*DELTA*((ALPHA*Z*MC)/U)^(1.00/(1.00-ALPHA)))/S)
        K          =       ((ALPHA*Z*MC)/U)^(1.00/(1.00-ALPHA))*H
        W          =       (1.00-ALPHA)*MC*Z*(K/H)^ALPHA
        YSS        =       Z*H^(1.00-ALPHA)*K^ALPHA
        

        PAR = [THETA; BETTA; GAMMA; DELTA; ALPHA; XI; RHO; PHI_R; PHI_PPI; PHI_Y; PPS; YSS; RRS; SIGMAA; SIGMAMP; SIGMAETA; RHOETA]

    # Solution functions (no adjustment needed)    
        eta     =   eval_ShockVAR(PAR)
        PAR_SS  =   eval_PAR_SS(PAR)
        SS      =   eval_SS(PAR_SS)
        deriv   =   eval_deriv(PAR_SS, SS)
        
        # @btime sol_mat = solve_model(model, deriv, eta)
        sol_mat = solve_model(model, deriv, eta)

# Store gx and hx

        gx=sol_mat[1];
        hx=sol_mat[2];

###################################################################################
## STEP N°2 = Define the IRFs function for the three shocks in the question      ##
###################################################################################

#IRF FOR TFP SHOCK
    function TFP(gx, hx, T)
            
        PVt = zeros(T, 5)
        NPVt = zeros(T, 13) 

        PVt[1, :] = [0, 1, 0, 0, 0] 

        for i in 2:T  
            NPVt[i - 1, :] = gx * PVt[i - 1, :]
            PVt[i, :] = hx * PVt[i - 1, :]
        end

        return NPVt, PVt
    end


#IRF FOR  MONETARY SHOCK
    function Monetary(gx, hx, T)
            
        PVt = zeros(T, 5) 
        NPVt = zeros(T, 13) 

        PVt[1, :] = [0, 0, 1, 0, 0] 

        for i in 2:T  
            NPVt[i - 1, :] = gx * PVt[i - 1, :]
            PVt[i, :] = hx * PVt[i - 1, :]
        end

        return NPVt, PVt
    end


#IRF FOR COST PUSH SHOCK
    function CostPush(gx, hx, T)
            
        PVt = zeros(T, 5) 
        NPVt = zeros(T, 13) 

        PVt[1, :] = [0, 0, 0, 0, 1] 

        for i in 2:T  
            NPVt[i - 1, :] = gx * PVt[i - 1, :]
            PVt[i, :] = hx * PVt[i - 1, :]
        end

        return NPVt, PVt
    end

########################################################################################
## STEP N°3 = Plot the IRF for general and relevant variables (Simulation for 3 cases)##
########################################################################################

#IRF FOR TFP SHOCK
    T = 50
    NPVt, PVt = TFP(gx, hx, T) 
    using Plots

    ## RELEVANT VARIABLES ##
    plot(1:T, NPVt[:, 1], label="consumption", color=:blue) 
    plot!(1:T, NPVt[:, 3], label="Inflation", color=:red)
    plot!(1:T, NPVt[:, 5], label="Policy rate", color=:green)
    plot!(1:T, NPVt[:, 8], label="Marginal cost of firms", color=:black)
    plot!(1:T, PVt[:, 1], label="Price dispersion", color=:orange)
    plot!(1:T, PVt[:, 2], label="TFP", color=:purple)
    plot!(1:T, PVt[:, 4], label="Capital", color=:gold)
    savefig("Figure 1 - Relevant Variables TFP.png")

    ##NOMINAL VARIABLES ## 
    plot(1:T, NPVt[:, 3], label="Inflation", color=:red)
    plot!(1:T, NPVt[:, 5], label="Policy rate", color=:green)
    plot!(1:T, NPVt[:, 8], label="Marginal cost of firms", color=:black)
    plot!(1:T, PVt[:, 1], label="Price dispersion", color=:orange)
    savefig("Figure 2 - Nominal Variables TFP.png")

#FIGURES INCLUDED IN THE ANNEX A ##########
    ##NON - PREDETERMINED VARIABLES##
    plot(1:T, NPVt[:, 1], label="consumption", color=:blue) 
    plot!(1:T, NPVt[:, 2], label="demand of hours of work")
    plot!(1:T, NPVt[:, 3], label="Inflation", color=:red)
    plot!(1:T, NPVt[:, 4], label="Price of asset")
    plot!(1:T, NPVt[:, 5], label="Policy rate", color=:green)
    plot!(1:T, NPVt[:, 6], label="Real wage", color=:black)
    plot!(1:T, NPVt[:, 7], label="Lagrangian multiplier hh")
    plot!(1:T, NPVt[:, 8], label="Marginal cost of firms", color=:black)
    plot!(1:T, NPVt[:, 9], label="Reoptimize Price level")
    plot!(1:T, NPVt[:, 10], label="price optimal setting 1")
    plot!(1:T, NPVt[:, 11], label="price optimal setting 2")
    plot!(1:T, NPVt[:, 12], label="Investment in Capital stock", color=:grey)
    plot!(1:T, NPVt[:, 13], label="Rental rate of capital")
    savefig("Figure 1A.png")

    ## PREDETERMINED VARIABLES ##
    plot(1:T, PVt[:, 1], label="Price dispersion", color=:orange)
    plot!(1:T, PVt[:, 2], label="TFP", color=:violet)
    plot!(1:T, PVt[:, 3], label="Target Policy rate")
    plot!(1:T, PVt[:, 4], label="Capital", color=:yellow)
    plot!(1:T, PVt[:, 5], label="Elasticty of demand")
    savefig("Figure 2A.png")

#IRF FOR  MONETARY SHOCK 
    T = 50
    NPVt, PVt = Monetary(gx, hx, T) 

    ###### RELEVANT VARIABLES ##########

    plot(1:T, NPVt[:, 1], label="consumption", color=:blue,) 
    plot!(1:T, NPVt[:, 3], label="Inflation", color=:red)
    plot!(1:T, NPVt[:, 4], label="Price of asset")
    plot!(1:T, NPVt[:, 5], label="Policy rate",color=:green)
    plot!(1:T, NPVt[:, 12], label="Investment in Capital stock", color=:grey)
    plot!(1:T, PVt[:, 1], label="Price dispersion", color=:orange)
    plot!(1:T, PVt[:, 3], label="Target Policy rate",color=:violet)
    plot!(1:T, PVt[:, 4], label="Capital",color=:yellow)
    savefig("Figure 3 - Relevant Variables Monetary Shock.png")

    ## NOMINAL VARIABLES ##
    plot(1:T, PVt[:, 3], label="Inflation", color=:red)
    plot!(1:T, PVt[:, 5], label="Policy rate",color=:green)
    plot!(1:T, NPVt[:, 1], label="Price dispersion", color=:orange)
    plot!(1:T, PVt[:, 4], label="Price of asset")
    savefig("Figure 4 - Nominal Variables Monetary Shock.png")

#FIGURES INCLUDED IN THE ANNEX A
    ## NON PREDETERMINED VARIABLES##
    plot(1:T, NPVt[:, 1], label="consumption", color=:blue) 
    plot!(1:T, NPVt[:, 2], label="demand of hours of work")
    plot!(1:T, NPVt[:, 3], label="Inflation", color=:red)
    plot!(1:T, NPVt[:, 4], label="Price of asset")
    plot!(1:T, NPVt[:, 5], label="Policy rate", color=:green)
    plot!(1:T, NPVt[:, 6], label="Real wage")
    plot!(1:T, NPVt[:, 7], label="Lagrangian multiplier hh")
    plot!(1:T, NPVt[:, 8], label="Marginal cost of firms", color=:black)
    plot!(1:T, NPVt[:, 9], label="Reoptimize Price level")
    plot!(1:T, NPVt[:, 10], label="price optimal setting 1")
    plot!(1:T, NPVt[:, 11], label="price optimal setting 2")
    plot!(1:T, NPVt[:, 12], label="Investment in Capital stock", color=:grey)
    plot!(1:T, NPVt[:, 13], label="Rental rate (price of capital)")
    savefig("Figure 3A.png")

    ## PREDETERMINED VARIABLES ##
    plot(1:T, PVt[:, 1], label="Price dispersion", color=:orange)
    plot!(1:T, PVt[:, 2], label="TFP")
    plot!(1:T, PVt[:, 3], label="Target Policy rate", color=:violet)
    plot!(1:T, PVt[:, 4], label="Capital", color=:yellow)
    plot!(1:T, PVt[:, 5], label="cost push(elasticity of demand w.r.t. prices)")
    savefig("Figure 4A.png")


#IRF FOR COST PUSH SHOCK
    T = 50
    NPVt, PVt = CostPush(gx, hx, T) 

    using Plots

    ## RELEVANT VARIABLES ## 
    plot(1:T, NPVt[:, 1], label="consumption", color=:blue)
    plot!(1:T, NPVt[:, 3], label="Inflation", color=:red)
    plot!(1:T, NPVt[:, 5], label="Policy rate",color=:green)
    plot!(1:T, NPVt[:, 8], label="Marginal cost of firms",color=:black)
    plot!(1:T, NPVt[:, 12], label="Investment in Capital stock", color=:grey)
    plot!(1:T, PVt[:, 1], label="Price dispersion", color=:orange)
    plot!(1:T, PVt[:, 4], label="Capital",color=:yellow)
    plot!(1:T, PVt[:, 5], label="cost push(elasticty of demand wrt prices)",color=:violet)
    savefig("Figure 7 - Relevant Variables Cost Push - Shock.png")

    ## NOMINAL VARIABLES ## 
    plot(1:T, NPVt[:, 3], label="Inflation", color=:red)
    plot!(1:T, NPVt[:, 5], label="Policy rate", color=:green)
    plot!(1:T, PVt[:, 1], label="Price dispersion", color=:orange)
    savefig("Figure 6 - Nominal Variables Cost Push - Shock.png")

    ###### FIGURES INCLUDED IN THE ANNEX A ##########

    ## UNPREDETERMINED VARIABLES ##
    plot(1:T, NPVt[:, 1], label="consumption", color=:blue)
    plot!(1:T, NPVt[:, 2], label="demand of hours of work")
    plot!(1:T, NPVt[:, 3], label="Inflation", color=:red)
    plot!(1:T, NPVt[:, 4], label="Price of asset")
    plot!(1:T, NPVt[:, 5], label="Policy rate", color=:green)
    plot!(1:T, NPVt[:, 6], label="Real wage")
    plot!(1:T, NPVt[:, 7], label="Lagrangian multiplier hh")
    plot!(1:T, NPVt[:, 8], label="Marginal cost of firms", color=:black)
    plot!(1:T, NPVt[:, 9], label="Reoptimize Price level")
    plot!(1:T, NPVt[:, 10], label="price optimal setting 1")
    plot!(1:T, NPVt[:, 11], label="price optimal setting 2")
    plot!(1:T, NPVt[:, 12], label="Investment in Capital stock", color=:grey)
    plot!(1:T, NPVt[:, 13], label="Rental rate (price of capital)")
    savefig("Figure 5A.png")

    ## PREDETERMINED VARIABLES ##
    plot(1:T, PVt[:, 1], label="Price dispersion", color=:orange)
    plot!(1:T, PVt[:, 2], label="TFP")
    plot!(1:T, PVt[:, 3], label="Target Policy rate")
    plot!(1:T, PVt[:, 4], label="Capital", color=:yellow)
    plot!(1:T, PVt[:, 5], label="cost push(elasticity of demand w.r.t. prices)", color=:violet)
    savefig("Figure 6A.png")







