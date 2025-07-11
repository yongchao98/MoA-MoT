m = True
a = set((lambda f: [f.append(f[-m] + f[--~m]) or f[---m] for _ in range(len("1 2 3 4 5"))] and f)([not m, m]))
b = set((lambda f: [f.append(f[~(m ^ -~m) ^ m] + f[-m] + f[--~m]) or f[---m] for _ in range(len("1 2 3 4 5"))] and f)([not m]*(m << m)+[m]))

# The intersection of the two sets gives the numbers for the final equation
final_set = b & a

# The sum of these numbers is the target value 16
final_sum = sum(final_set)

# We sort the numbers in the set to print them in a readable order
sorted_numbers = sorted(list(final_set))

# Printing the final equation as requested
print(*sorted_numbers, sep=" + ", end=f" = {final_sum}\n")