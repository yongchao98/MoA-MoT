m = True
a = set((lambda f: [f.append(f[-m] + f[--~m]) or f[---m] for _ in range(len("1 2 3 4 5"))] and f)([not m, m]))
b = set((lambda f: [f.append(f[~(m ^ -~m) ^ m] + f[-m] + f[--~m]) or f[---m] for _ in range(len("1 2 3 4 5"))] and f)([not m]*(m << m)+[m]))

# The final operation is sum(b & a)
# The resulting set from the intersection is b & a
result_set = b & a
# The final result is the sum of the elements in the intersection set
final_result = sum(result_set)

# To print the equation, we sort the numbers and join them with ' + '
equation_parts = sorted(list(result_set))
print(f"{' + '.join(map(str, equation_parts))} = {final_result}")