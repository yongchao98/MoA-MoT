# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over possible values of B
for B in numbers:
    A = B + 9
    J = B + 3
    C = 52 - B
    K = 57 - B
    L = B + 17
    H = 12 - B
    E = (52 - B) / 3
    D = E - 13
    I = 18 - E
    G = 3.6 * J
    
    # Check if all values are integers and in the list of given numbers
    if all(x in numbers for x in [A, B, C, D, E, F, G, H, I, J, K, L, M]):
        # If all conditions are satisfied, print the result
        print(f"<<<{[A, B, C, D, E, F, G, H, I, J, K, L, M]}>>>")
        break