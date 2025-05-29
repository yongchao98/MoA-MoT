# Possible values for the letters
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over possible values for A
for A in values:
    D = A + 43
    F = 43 - A
    if F % 1.5 != 0:
        continue
    C = F / 1.5
    if C % 1.6 != 0:
        continue
    L = C / 1.6
    if L % 3.0 != 0:
        continue
    J = L / 3.0
    if F - J != 31:
        continue
    if F % 3.6 != 0:
        continue
    G = F / 3.6
    M = 4.5 * G
    E = 12 - G
    K = 33 - J
    
    # Check if all calculated values are in the list and distinct
    calculated_values = [A, C, D, E, F, G, J, K, L, M]
    if all(val in values for val in calculated_values) and len(set(calculated_values)) == len(calculated_values):
        # Find remaining values for B, H, I
        remaining_values = set(values) - set(calculated_values)
        if len(remaining_values) == 3:
            B, H, I = remaining_values
            if H > L and M > A:
                # Print the result in alphabetical order
                result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
                print(result)
                break