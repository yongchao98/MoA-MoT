import itertools

# List of numbers to assign
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the numbers
for perm in itertools.permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all conditions
    if (C - G == -36 and
        C - F == -17 and
        B + F == 68 and
        E + H == 137 and
        A + K == 375 and
        C + H == 124 and
        A > F and
        J > C and
        F - G == -19 and
        A > C and
        C + G == 42):
        
        # Output the result in alphabetical order
        result = [A, B, C, D, E, F, G, H, I, J, K]
        print(f"<<<{result}>>>")
        break