########## R function: updateTwoLevGaussPenFrag ##########

# For updating the message graph and natural parameter vectors of the
# two-level linear mixed model Gaussian penalization fragment.

# Last changed: 18 MAY 2020

updateGaussPen2levFrag <- function(mu.beta,Sigma.beta,etaIN)
{
   # The naming order for etaIN is:

   # 1. p(beta,u|Sigma)->(beta,u)
   # 2. (beta,u)->p(beta,u|Sigma)
   # 3. p(beta,u|Sigma)->Sigma
   # 4. Sigma->p(beta,u|Sigma)

   # The naming order for etaOUT is:

   # 1. "p(beta,u|Sigma)->(beta,u)",
   # 2. "p(beta,u|Sigma)->Sigma"

   # Obtain dimension variables:

   p <- length(mu.beta)
   lenVal <- length(etaIN[[3]]) - 1
   q <- (sqrt(8*lenVal + 1) - 1)/2

   m <- (length(etaIN[[1]]) - p - 0.5*p*(p + 1))/(q + 0.5*q*(q + 1) + p*q)

   # Set required duplication matrices and their Moore-Penrose
   # inverses:

   Dp <- duplication.matrix(p)  ;   DpPlus <- solve(crossprod(Dp),t(Dp))
   Dq <- duplication.matrix(q)  ;   DqPlus <- solve(crossprod(Dq),t(Dq))

   # Allocate natural parameter lists:

   etaSUM <- vector("list",2)  ;  etaOUT <- vector("list", 2)

   # Obtain the sums of natural parameters along each edge:

   etaSUM$"p(beta,u|Sigma)<->(beta,u)" <- etaIN[[1]] + etaIN[[2]]
   etaSUM$"p(beta,u|Sigma)<->Sigma" <- etaIN[[3]] + etaIN[[4]]

   # Update expectations required for the factor to stochastic
   # node message natural parameters:

   omega23 <- etaSUM$"p(beta,u|Sigma)<->Sigma"[1]
   omega24 <- etaSUM$"p(beta,u|Sigma)<->Sigma"[-1]

   M.q.inv.Sigma <- (omega23 + 0.5*(q+1))*solve(vecInverse(crossprod(DqPlus,omega24)))

   # Convert etaSUM$"p(y|beta,u,Sigma<->(beta,u)" to natural parameters
   # of interest:

   eta_in_nu <- etaIN[1:2]
   
   S5 <- two_level_nat_to_comm_parms(p, q, m, eta_in_nu)

   # Obtain summation over groups required for covariance matrix parameter natural parameter:

   omega25 <- rep(0,0.5*q*(q+1))
   for (i in 1:m)
      omega25 <- omega25 - 0.5*crossprod(Dq,vec(tcrossprod(S5$E_q_nu_2[[i]]) + S5$Cov_q_nu_2[[i]]))

   # Obtain output lists:

   etaOUT1 <- c(solve(Sigma.beta,mu.beta),-0.5*crossprod(Dp,vec(solve(Sigma.beta))))
   for (i in 1:m)
   {
      currVec <-  c(rep(0,q),-0.5*crossprod(Dq,vec(M.q.inv.Sigma)),rep(0,p*q))
      etaOUT1 <- c(etaOUT1,currVec)
   }

   etaOUT[[1]] <- etaOUT1

   etaOUT[[2]] <- c((-0.5*m),omega25)

   return(etaOUT)
}

############ End of updateTwoLevGaussPenFrag ############
