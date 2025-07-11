m = True

# With index -1: new = f[-1] + f[-1] = 2 * f[-1]
# A=12 (-), f starts [-1, 1]. seq: 1, 2, 4, 8, 16, 32, ...
a_list = [-1, 1]
for _ in range(len("1 2 3 4 5")):
    a_list.append(a_list[-1] + a_list[-1])
a = set(a_list) # a = {-1, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512}

# With index -1: idx = -1, f[1]=-1. new = f[-1] + f[-1] + f[-1] = 3*f[-1]
# B=2 (>>), f starts [-1,-1,1].
b_list = [-1, -1, 1]
for _ in range(len("1 2 3 4 5")):
    # idx from B=>> is -1. index in f.append is m -> -1.
    b_list.append(b_list[-1] + b_list[-1] + b_list[-1]) # using index -1 for both expressions
b = set(b_list) # b = {-1, 1, 3, 9, 27, ...}

# F=6(^), E=11(sum). sum(a^b) = sum({2,3,4,8,9,16,27...}) != 16.