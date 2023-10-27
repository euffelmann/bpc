h2o_to_h2l<-function(K,P = 0.5,h2o = 1){
  ## Transform observed scale h^2 to liability scale h^2
  ## Based on Eq 23 of Lee et al 2011 AJHG
  t = -qnorm(K,0,1) ; z = dnorm(t)
  return(h2o*K*K*(1-K)*(1-K)/{z*z*P*(1-P)})
}
