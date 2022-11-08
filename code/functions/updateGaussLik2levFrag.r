########## R function: updateTwoLevGaussLikFrag ##########

# For updating the message graph and natural parameter vectors of the
# two-level linear mixed model Gaussian likelihood fragment.

# Last changed: 18 MAY 2020

updateGaussLik2levFrag <- function(yList,XList,Zlist,etaIN)
{
   # The naming order for etaIN is:

   # 1. p(y|beta,u,sigsq)->(beta,u)
   # 2. (beta,u)->p(y|beta,u,sigsq)
   # 3. p(y|beta,u,sigsq)->sigsq
   # 4. sigsq->p(y|beta,u,sigsq)

   # The naming order for etaOUT is:

   # 1. "p(y|beta,u,sigsq)->(beta,u)",
   # 2. "p(y|beta,u,sigsq)->sigsq"

   # Obtain dimension variables:

   m <- length(yList)
   p <- ncol(XList[[1]])   ;   q <- ncol(Zlist[[1]])

   nVec <- as.numeric(unlist(lapply(y,length)))

   # Set required duplication matrices and their Moore-Penrose
   # inverses:

   Dp <- duplication.matrix(p)  ;   DpPlus <- solve(crossprod(Dp),t(Dp))
   Dq <- duplication.matrix(q)  ;   DqPlus <- solve(crossprod(Dq),t(Dq))

   # Allocate natural parameter lists:

   etaSUM <- vector("list",2) ; etaOUT <- vector("list",2)

   # Obtain the sums of natural parameters along each edge:

   etaSUM$"p(y|beta,u,sigsq)<->(beta,u)" <- etaIN[[1]] + etaIN[[2]]
   etaSUM$"p(y|beta,u,sigsq)<->sigsq" <- etaIN[[3]] + etaIN[[4]]

   # Update expectations required for the factor to stochastic
   # node message natural parameters:

   mu.q.recip.sigsq <- ((etaSUM$"p(y|beta,u,sigsq)<->sigsq"[1] + 1)
                           /(etaSUM$"p(y|beta,u,sigsq)<->sigsq"[2]))

   # Convert etaSUM$"p(y|beta,u,sigsq)<->(beta,u)" to natural parameters
   # of interest:
   
   eta_in_nu <- etaIN[1:2]

   S4 <- two_level_nat_to_comm_parms(p, q, m, eta_in_nu)

   mu.q.beta <- S4$"E_q_nu_1"  ;  Sigma.q.beta <- S4$"Cov_q_nu_1"
   mu.q.u <- S4$"E_q_nu_2"  ;  Sigma.q.u <- S4$"Cov_q_nu_2"
   E.q.CROSS.beta.u <- S4$"Cov_q_nu_12"

   # Obtain summation over groups required for variance parameter natural parameter:

   omega20 <- rep(0,p)   ;  omega21 <- rep(0,(0.5*p*(p+1))) ; omega22 <- 0
   for (i in 1:m)
   {
      omega20 <- omega20 + crossprod(XList[[i]],yList[[i]])
      omega21 <- omega21 - 0.5*crossprod(Dp,vec(crossprod(XList[[i]])))
      residVec <- (yList[[i]] - as.vector(crossprod(t(XList[[i]]),mu.q.beta))
                            - as.vector(crossprod(t(Zlist[[i]]),mu.q.u[[i]])))
      omega22 <- omega22 - 0.5*sum(residVec^2)
      omega22 <- omega22 - 0.5*sum(diag(crossprod(Sigma.q.beta,crossprod(XList[[i]]))))
      omega22 <- omega22 - 0.5*sum(diag(crossprod(Sigma.q.u[[i]],crossprod(Zlist[[i]]))))
      omega22 <- omega22 - sum(diag(crossprod(crossprod(XList[[i]],Zlist[[i]]),E.q.CROSS.beta.u[[i]])))
   }

   etaOUT1 <- mu.q.recip.sigsq*c(omega20,omega21)
   for (i in 1:m)
   {
      currVec <-  c(crossprod(Zlist[[i]],yList[[i]]),
                    -0.5*crossprod(Dq,vec(crossprod(Zlist[[i]]))),
                    -vec(crossprod(XList[[i]],Zlist[[i]])))
      etaOUT1 <- c(etaOUT1,mu.q.recip.sigsq*currVec)
   }

   etaOUT[[1]] <-  etaOUT1
   etaOUT[[2]] <-  c(-0.5*sum(nVec),omega22)

   return(etaOUT)
}

############ End of updateTwoLevGaussLikFrag ############
