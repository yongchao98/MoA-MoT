x = 11**8
while 1:
    x = (x*(2**25-x))>>25
    print '1' if x>>(25-1) else '0',