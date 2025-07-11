k,a,b,a1,b1=2,4,1,12,4
while 1:
 p,q,k=k*k,2*k+1,k+1
 a,b,a1,b1 = a1, b1, p*a+q*a1, p*b+q*b1
 d,d1=a/b, a1/b1
 while d==d1:
  print(int(d))
  a,b=(a%b)*2,(b)
  a1,b1=(a1%b1)*2,(b1)
  d,d1=a/b,a1/b1