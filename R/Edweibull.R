Edweibull <-
function(q, beta, eps=0.0001, nmax=1000, zero=FALSE)
{
if(beta==1 & !zero)
1/(1-q)
else if (beta==1 & zero)
q/(1-q)
else
{
xmax<-min(2*qdweibull(1-eps, q, beta, zero),nmax)
x<-1:xmax
sum(ddweibull(x, q, beta, zero)*x)
}
}

E2dweibull <-
function(q, beta, eps=0.0001, nmax=1000, zero=FALSE)
{
if(beta==1 & !zero)
(1+q)/(1-q)^2
if(beta==1 & zero)
q*(1+q)/(1-q)^2
else
{
xmax<-min(2*qdweibull(1-eps, q, beta, zero),nmax)
x<-1:xmax
sum(ddweibull(x, q, beta, zero)*x^2)
}
}

Vdweibull <-
function(q, beta, eps=0.0001, nmax=1000, zero=FALSE)
{
if(beta==1 & !zero)
q/(1-q)^2
else
{
xmax<-min(2*qdweibull(1-eps, q, beta, zero),nmax)
E2dweibull(q, beta, eps, xmax, zero)-Edweibull(q, beta, eps, xmax, zero)^2
}
}

ERdweibull <-
function(q, beta, eps=0.0001, nmax=1000)
{
if(beta==1)
(1-q)/q*log(1/(1-q))
else
{
xmax<-min(2*qdweibull(1-eps, q, beta),nmax)
x<-1:xmax
sum((q^(x-1)^beta-q^x^beta)*1/x)
}
}

