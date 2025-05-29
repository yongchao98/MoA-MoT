# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Known values
A = 5
F = 10
L = 16
E = 9

# Iterate over possible values for J
for J in numbers:
    I = 4 * J
    H = I + 4
    M = 22 - 4 * J
    D = 6 * J + 6
    B = D - 14
    
    # Check if all values are unique and within the given numbers
    if (I in numbers and H in numbers and M in numbers and D in numbers and B in numbers and
        len(set([A, B, C, D, E, F, G, H, I, J, K, L, M])) == 13):
        # Assign remaining values
        remaining_numbers = set(numbers) - {A, B, D, E, F, H, I, J, L, M}
        C, G, K = remaining_numbers
        
        # Print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break