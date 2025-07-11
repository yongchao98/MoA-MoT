x = 1
K = 0x39048a39ec18a4a6
# The exec forms a loop to generate the sequence
exec("y = x & 1; print(y); x = (x >> 1) ^ ([0, K][y])")