import sys
n=int(sys.argv[1])
l=[1,2,2]
for i in range(2,n):l+=[1+i%2]*l[i]
print(l[:n])