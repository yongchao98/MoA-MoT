# List of possible values
values = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Find possible values for K, D, and C
for K in values:
    D = 3 * K
    C = 2 * K
    if D in values and C in values:
        print(f"K: {K}, D: {D}, C: {C}")