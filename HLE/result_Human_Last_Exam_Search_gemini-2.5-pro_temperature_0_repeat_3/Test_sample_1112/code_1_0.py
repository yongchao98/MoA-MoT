def Kolakoski():
    x = y = -1
    while True:
        yield [2,1][x&1]
        f = y &~ (y+1)
        x ^= f
        y = (y+1) | (f & (x>>1))