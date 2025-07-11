a=[1,2];i=s=0
while 1:
 s^=1;print(end=['1','0'][s]);i+=1
 if i==len(a):a+=a[-1]*[a[i-2]]