# Given numbers
numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Initialize a dictionary to store the values of each letter
values = {}

# Solve for G and E
for G in numbers:
    E = G - 27
    if E in numbers:
        # Solve for C
        C = E - 5
        if C in numbers:
            # Solve for H
            H = 4 * C
            if H in numbers:
                # Verify E + H = 25
                if E + H == 25:
                    # Solve for A
                    A = 137 - H
                    if A in numbers:
                        # Solve for D and F
                        for D in numbers:
                            F = D - 79
                            if F in numbers:
                                # Check if all values are unique
                                if len(set([A, C, D, E, F, G, H])) == 7:
                                    # Find the remaining number for B and I
                                    remaining_numbers = set(numbers) - {A, C, D, E, F, G, H}
                                    B, I = remaining_numbers
                                    # Assign values to the dictionary
                                    values['A'] = A
                                    values['B'] = B
                                    values['C'] = C
                                    values['D'] = D
                                    values['E'] = E
                                    values['F'] = F
                                    values['G'] = G
                                    values['H'] = H
                                    values['I'] = I
                                    break

# Output the values in alphabetical order
result = [values[letter] for letter in sorted(values.keys())]
print(f"<<<{result}>>>")