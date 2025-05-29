def bf(c):
    m=[0]*30000;p=i=0;o=''
    while i<len(c):
        if c[i]=='>':p+=1
        elif c[i]=='<':p-=1
        elif c[i]=='+':m[p]=(m[p]+1)%256
        elif c[i]=='-':m[p]=(m[p]-1)%256
        elif c[i]=='.':o+=chr(m[p])
        elif c[i]=='[' and not m[p]:
            n=1
            while n:i+=1;n+=(c[i]=='[')-1*(c[i]==']')
        elif c[i]==']' and m[p]:
            n=1
            while n:i-=1;n+=(c[i]==']')-1*(c[i]=='[')
            continue
        i+=1
    print(o,end='')

bf("[-]>[-]<>++++++++[<++++++++++++++>-]<++.-------------.+++++++++++++.+++.-------.<")