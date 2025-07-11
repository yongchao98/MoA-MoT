s='1'
i=0
while 1:
    t=s.count('1',i)
    s+='0'*t+bin(len(s)-i)[2:]+'1'
    i+=t+len(bin(len(s)-i-t)[2:])