from itertools import permutations

# The numbers to assign
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all the conditions
    if (B + E == 170 and
        H - C == -46 and
        K == 3.0 * G and
        B + I == 230 and
        C - H == 46 and
        D == 3.0 * H and
        C - D == -104 and
        H - J == 36 and
        A == 2.4 * E and
        I - E == 60 and
        C - K == 112):
        
        # If all conditions are satisfied, print the result
        result = [A, B, C, D, E, F, G, H, I, J, K]
        print(f"<<<{result}>>>")
        break