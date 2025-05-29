# Given numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Initialize a dictionary to store the values of each letter
values = {}

# Solve the equations step by step
# From equation 3: C = B - 44
# From equation 9: D = 84 - C
# From equation 7: E = 57 - D
# From equation 5: I = 130 - E
# From equation 6: J = I + 29
# From equation 1: C = 186 - J
# From equation 8: K = B - 64
# From equation 10: G = B - 53
# From equation 4: H = 234 - E
# From equation 2: F = H / 3

# Start with B, since it appears in multiple equations
for B in numbers:
    C = B - 44
    D = 84 - C
    E = 57 - D
    I = 130 - E
    J = I + 29
    K = B - 64
    G = B - 53
    H = 234 - E
    F = H / 3

    # Check if all values are unique and in the given numbers
    if len(set([B, C, D, E, F, G, H, I, J, K])) == 10 and all(x in numbers for x in [B, C, D, E, F, G, H, I, J, K]):
        values = {'A': None, 'B': B, 'C': C, 'D': D, 'E': E, 'F': F, 'G': G, 'H': H, 'I': I, 'J': J, 'K': K}
        break

# Find the remaining number for A
remaining_number = list(set(numbers) - set(values.values()))[0]
values['A'] = remaining_number

# Extract the solution in alphabetical order
result = [values[letter] for letter in sorted(values.keys())]

# Print the result
print(f"<<<{result}>>>")