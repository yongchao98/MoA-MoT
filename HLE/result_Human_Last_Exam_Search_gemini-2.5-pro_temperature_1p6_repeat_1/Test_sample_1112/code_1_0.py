a,b,c,d,e,f=1,0,1,1,3,3
while 1:
 if 4*a+b-c<e*c:
  print(e)
  e,f=10*(a-e*c),10*f
 else:
  a,b,c,d,e,f=a*d, (2*a+b)*f,c*f,d+1, (a*(7*d+2)+b*f)//(c*f),f+2