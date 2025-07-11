m = True
a = set((lambda f: [f.append(f[-m] + f[--~m]) or f[---m] for _ in range(len("1 2 3 4 5"))] and f)([not m, m]))
b = set((lambda f: [f.append(f[~(m ^ -~m) ^ m] + f[-m] + f[--~m]) or f[---m] for _ in range(len("1 2 3 4 5"))] and f)([not m]*(m << m)+[m]))

# The code above calculates two sets:
# a = {0, 1, 2, 3, 5, 8, 13, 21, 34, 55} (Fibonacci numbers)
# b = {0, 1, 2, 4, 7, 13, 24, 44, 81, 149} (Tribonacci numbers)

# The final operation is sum(b & a).
# The intersection of the two sets is {0, 1, 2, 13}.
# The sum of these elements is 0 + 1 + 2 + 13 = 16.

final_set = b & a
# The final part of the code prints the equation as requested.
print(f'{" + ".join(map(str, sorted(list(final_set))))} = {sum(final_set)}')
