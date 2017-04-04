
library("EpiModel")

####################################################################
##################### Plot Basic SI Model ##########################
####################################################################
param <- param.dcm(inf.prob = 0.5, act.rate = 0.25)
init <- init.dcm(s.num = 1000, i.num = 1)
control <- control.dcm(type = "SI", nsteps = 150)

SImodel <- dcm(param, init, control)

plot(SImodel, main="Basic SI Model")



####################################################################
###################SIR NO DEMOGRAPHICS##############################
####################################################################

# Write Model For SIR No Demographic Example: Influenza at a Boarding School 
SIRND <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    # Population size
    num <- s.num + i.num + r.num
    
    # Effective contact rate and FOI from a rearrangement of Beta * c * D
    ce <- R0 / i.dur #Beta
    lambda <- ce * i.num/num #Beta*I
    
    dS <- -lambda*s.num #-Beta*S*I
    dI <-  lambda*s.num - (1/i.dur)*i.num 
    dR <- (1/i.dur)*i.num
    
    # Compartments and flows are part of the derivative vector
    # Other calculations to be output are outside the vector, but within the containing list
    list(c(dS, dI, dR, 
           si.flow = lambda*s.num,
           ir.flow = (1/i.dur)*i.num,
         num = num,
        i.prev = i.num/num,
         s.prev = s.num/num))
  })
}

#Plot SIR No Demographics Example Influenza in a Boarding School 

param <- param.dcm(R0 = 3.65, i.dur = 2.2)
init <- init.dcm(s.num = 760, i.num = 3, r.num = 0, si.flow = 0, ir.flow = 0, num = 763, i.prev = 3/763, s.prev = 760/763)
control <- control.dcm(nsteps = 15, dt = 1, new.mod = SIRND)
mod <- dcm(param, init, control)

plot(mod, y = c("s.num", "i.num"), main = "SIR No Demographics: Influenza in Boarding School", leg = "full")




##################################################
#########SIR WITH DEMOGRAPHICS EXAMPLE############
##################################################

#Write SIR With Demographic Model

SIRWD <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    # Population size
    num <- s.num + i.num + r.num
    
    # Effective contact rate and FOI from a rearrangement of Beta * c * D
    ce <- R0 / i.dur #Beta
    lambda <- ce * i.num/num #Beta * I 
    s.prev <- s.num/num
    i.prev <- i.num/num
    r.prev <- r.num/num
    
    
    dS <- b.rate - lambda*s.num - ds.rate*s.num
    dI <-  lambda*s.num - (1/i.dur)*i.num - di.rate*i.num
    dR <- (1/i.dur)*i.num - dr.rate*r.num
    
    # Compartments and flows are part of the derivative vector
    # Other calculations to be output are outside the vector, but within the containing list
    list(c(dS, dI, dR, 
           si.flow = lambda*s.num,
           ir.flow = (1/i.dur)*i.num),
        num = num,
        i.prev = i.num/num,
        s.prev = s.num/num)
  })
}

# Input Parameters and Plot
param <- param.dcm(R0 = 5, i.dur = 1, b.rate = 1/70, ds.rate = 1/70, di.rate = 1/70, dr.rate = 1/70)
init <- init.dcm(s.num = 10000, i.num = 25, r.num = 0, si.flow = 0, ir.flow = 0)
control <- control.dcm(nsteps = 1000, dt = 1, new.mod = SIRWD)
mod <- dcm(param, init, control)

plot(mod, y = "i.prev", ylim = c(0, 0.5), main = "SIR With Demographics", leg = "full")

#, num = 100000, i.prev = 2.5e-4, s.prev = 0.1)
##########SIR WITH DEMOGRAPHY EPIMODEL EXAMPLE###########################

param <- param.dcm(inf.prob = 0.2, act.rate = 1, rec.rate = 1/20,
                   b.rate = 1/70, ds.rate = 1/70, di.rate = 1/70, dr.rate = 1/70)
init <- init.dcm(s.num = 1000, i.num = 2, r.num = 0)
control <- control.dcm(type = "SIR", nsteps = 500, dt = 0.5)
mod <- dcm(param, init, control)

plot(mod, main = "SIR With Demography: Epimodel Example")


####################################################################
########################### SEIR MODEL##############################
####################################################################

#######Write SEIR Model#########
SEIR <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    # Population size
    num <- s.num + e.num + i.num + r.num
    
    # Effective contact rate and FOI from a rearrangement of Beta * c * D
    ce <- R0 / i.dur #beta
    lambda <- ce * i.num/num #beta*I
    
    dS <- -lambda*s.num #-B*S*I
    dE <- lambda*s.num - (1/e.dur)*e.num 
    dI <- (1/e.dur)*e.num - (1 - cfr)*(1/i.dur)*i.num - cfr*(1/i.dur)*i.num #cfr = P(disease inducde mortality)
    dR <- (1 - cfr)*(1/i.dur)*i.num
    
    # Compartments and flows are part of the derivative vector
    # Other calculations to be output are outside the vector, but within the containing list
    list(c(dS, dE, dI, dR, 
           se.flow = lambda * s.num,
           ei.flow = (1/e.dur) * e.num,
           ir.flow = (1 - cfr)*(1/i.dur) * i.num,
           d.flow = cfr*(1/i.dur)*i.num),
         num = num,
         i.prev = i.num / num,
         s.prev = s.num/num)
  })
}

####Input Parameters and Plot#####

param <- param.dcm(R0 = 1.9, e.dur = 10, i.dur = 14, cfr = 0)
init <- init.dcm(s.num = 1e6, e.num = 10, i.num = 0, r.num = 0,
                 se.flow = 0, ei.flow = 0, ir.flow = 0, d.flow = 0)
control <- control.dcm(nsteps = 750, dt = 1, new.mod = SEIR)
mod <- dcm(param, init, control)

plot(mod, y = c("i.prev", "s.prev"), 
     main = "SEIR Model", leg = "full")



print(mod_SIR_1g_cl)



