m = True
a = set((lambda f: [f.append(f[-m] + f[--~m]) or f[---m] for _ in range(len("1 2 3 4 5"))] and f)([not m, m]))
b = set((lambda f: [f.append(f[~(m | -~m) | m] + f[-m] + f[--~m]) or f[---m] for _ in range(len("1 2 3 4 5"))] and f)([not m]*(m << m)+[m]))

# The expression to be calculated is sum(b & a)
# We find the numbers in the intersection and format the output as an equation.
intersection_set = b & a
numbers_to_sum = sorted(list(intersection_set))
equation_str = " + ".join(map(str, numbers_to_sum))
result = sum(numbers_to_sum)
print(f"{equation_str} = {result}")