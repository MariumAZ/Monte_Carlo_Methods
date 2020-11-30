#Exercice 1 :
#1.réponse sur pdf joint.
#2.
 M<-3*sqrt(2*pi)
 #Definition de la densite instrumentale:
 g<-function(x,y){
   return(dnorm(x,2,1)*dunif(y,0,2))
 }

 #generateur associe:
   rgen<-function(n){
       return(cbind (rnorm(n,2,1),runif(n,0,2)))
     }

 
 #densite de f1 (f tilde):
 f1<-function(x,y){
   a<-exp(-0.5*((x-2)**2))
   b<-(cos(x)**2)+2*(sin(y)**2)*(cos(x)**4)
   d<-1+4*((y-1)**2)
   return ((a+b)/d)
   
 }
 
#3.
#On simule f suivant la méthode de rejet:
 rejet<- function(n, rgen, g, f1, M) {
   ans <- matrix(NA, n, 2) # Echantillon de sortie
   counter <- 0 # Nombre de tirages suivant g
   for (i in 1:n) {
     z <- rgen(1) ; counter <- counter + 1
     while (runif(1) > (f1(z[,1],z[,2])/(M * g(z[,1],z[,2])))) {
       z <- rgen(1) ; counter <- counter + 1
     }
     ans[i, ] <- z
   }
   return(list(S = ans, N = counter))
 }


#4.
 #etape1: une fonction metropolis:
 Metropolis <- function(n,  log.g, log.f1) {
   x <- matrix(0, n+1, 2)
   # Initialisation
   x[1,] <- rgen(1) 
   r <- log.g(x[,1],x[,2]) - log.f1(x[,1],x[,2])
   for (i in 1:n) {
     # Proposition
     x.new <- rgen(1)
     r.new <- log.g(x.new[,1],x.new[,2]) - log.f1(x.new[,1],x.new[,2])
     if (log(runif(1)) <= r - r.new) { # Acceptation rejet
       x[i+1, ] <- x.new ; r <- r.new
     } else {
       x[i+1, ] <- x[i, ]
     }
   }
   return(x)
 }
 # etape 2 :Simulation d'une chaine de Markov de loi stationnaire
 log.f1<-function(x,y) log(f1(x,y))
 log.g<-function(x,y) log(g(x,y))

 

 ##exercice2:
 
 #partie I :
 
 #1.a/
 
 # On introduit la fonction Mc estim qui prend pour en  argument un vecteur
 #V contenant les 
 #realisations de variables aleatoires iid, un niveau de confiance level 
 #et retourne un data.frame:
 MC.estim<-function(v,level=0.95){
   n<-length(v) ; value<-cumsum(v)/(1:n) 
   s2<-(cumsum(v^2)-(1:n)*(value)^2)/(0:(n-1)) ; s2<-s2/(1:n)
   s2<-c(0, s2[-1])
   b.IC<-qnorm(0.5*(level+1))*sqrt(s2)
   return(data.frame(value=value , var =s2 , binf = value - b.IC , bsup= value + b.IC))
 }

 #La  methode classique:
 n=10000
 t=2
 x=rweibull(10000,2,1) #où x = x1+x2+x3
 h<-function(x){return(x>=t)}
 mc1<-h(x)
 I1<-MC.estim(h(x))
 
 ##tracer l'évolution
 plot(1:n, I1$value,type="l",col="steelblue4",lwd=1,xlab="n",ylim=c(0.01,0.03),main=expression(paste(bold("Estimation de "),italic(I),bold(" : m?thode classique"))))
 lines(2:n, I1$binf[-1],col="gold")
 lines(2:n, I1$bsup[-1],col="gold")
 legend("topright",c("Estimation","Intervalle de confiance à  95%"),
        col=c("steelblue4","gold"),
        lty=c(1,1),lwd=c(2,2),inset=0.05,bg="gray95",box.lty=0)
 
 #1.b/
 
 ##estimateur stratifié:
 
 MC.estim.strat<-function(n,L,level=0.95){
   n.l<-n/L
   # Realisations des lois conditionnelles pour X
   u<-runif(n)
   x1<-qnorm((rep(0:(L-1), n.l)+u)/L,0,2)
   
   ### On regarde l'evolution à chaque fois que l'on
   ### ajoute une observation à  toutes les strates
   y<-h(x)
   y.div<-split(y,as.factor(rep(1:L, n.l)))
   # Evolution de la moyenne par strate
   mean.s<-sapply(y.div, cumsum)/(1:n.l)
   # evolution de lestimation
   value<-rowSums(mean.s)/L
   # Evolution des variances intra-strates
   s2<-sapply(y.div,function(x)cumsum(x^2))
   s2<-(s2-(1:n.l)*mean.s^2)/(0:(n.l-1))
   # evolution de la variance
   s2.tot<-rowSums(s2)/(L*L*(1:n.l))
   b.IC<-qnorm(0.5*(level+1))*sqrt(s2.tot)
   return(data.frame(value= value,var= s2.tot,binf= value-b.IC,bsup= value+b.IC))}

  ### Experience 1 : 2 strates
 L<-2; p9.2<-MC.estim.strat(n, L)
 p9.2[n/L, ]
 
  ###Experience 2 : 10 strates
 L<-10; p9.10<-MC.estim.strat(n, L)
 p9.10[n/L, ]
 ### Experience 3 : 1000 strates
 L<-100; p9.100<-MC.estim.strat(n, L)
 p9.100[n/L, ]
 
 
 #2.a/
 
 #D'abord on simule suivant la loi de X3 en utilisant 
 #la methode de la fonction inverse de F:
 u<-runif(1)
 inv<-function(u){
   if (u>=0 & u<0.25)
     return (4*u)
   if (u>=0.25 & u<0.75)
     return (1)
   if (u>=0.75 & u<=1)
     return (4*u-2)
 }

 #on genere un echantillon  suivant la loi de X3:
 
 rgen2<-function(n){
  
  ans<-c() #echantillon de sortie
  for (i in 1:n)
  { u1<-runif(1)
    ans[i]=inv(u1)
  }
  return (ans) 
    
}
#estimateur par la methode classique :
n<-100000
x1<-rexp(n,1)
x2<-rexp(n,1)
x3<-rgen2(n)
h<-function(a)(a >=1)
mc<-h(x1+x2+x3)
p1<-MC.estim(mc,level=0.95)
#traçage de l'évolution de l'estimateur par la méthode classique:
plot(1:n, p1$value, type = "l", col = "blue4", lwd = 2,
     xlab = "n", ylim = c(0.95, 1), main = expression(paste(
       bold("Estimation de "), italic(p),bold(" : méthode classique"))))
lines(2:n, p1$binf[-1], col = "indianred1")

lines(2:n, p1$bsup[-1], col = "indianred1")
legend("topright", c("Estimation", "Intervalle de confiance à 95%"),
       col = c("blue4", "indianred1"), lty = c(1, 1), lwd = c(2, 2),
       inset = 0.05, bg = "gray95", box.lty = 0)

#2.b/reponse sur pdf joint.

 

  