a,b,c,d,e,f=0,1,1,0,1,2
while 1:
 if 4*a+3*b+c-d-e-f>0:
  g=int((4*a+3*b+c-d-e-f)/(a+b+c))
  a,b,c,d,e,f=(10*a,10*b,10*c,10*d,10*e,10*f,g)
  print("%x"%g,end="")
 else:
  a,b,c,d,e,f=(2*a+b,4*b+c,8*c+d,16*d+e,32*e+f,f+1)