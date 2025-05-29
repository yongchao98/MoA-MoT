# Given numbers
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Initialize a dictionary to store the values of the letters
values = {}

# Solve the constraints step-by-step
# From equation 9: B = 3.0F, so B must be a multiple of 3
# From equation 4: E + F = 42
# From equation 2: F + I = 83
# From equation 10: C + F = 23
# From equation 1: I = 4.0C

# Let's start by finding F, B, C, I, and E
for F in numbers:
    B = 3 * F
    if B not in numbers:
        continue
    for C in numbers:
        if C + F != 23:
            continue
        I = 4 * C
        if I not in numbers:
            continue
        if F + I != 83:
            continue
        E = 42 - F
        if E not in numbers:
            continue
        # Now we have F, B, C, I, E
        # Let's find K using equation 3: K - E = 82
        K = E + 82
        if K not in numbers:
            continue
        # Now we have F, B, C, I, E, K
        # Let's find A and H using equation 5: A + H = 225
        for A in numbers:
            H = 225 - A
            if H not in numbers:
                continue
            # Now we have A, H
            # Check equation 8: J = 3.0H
            J = 3 * H
            if J not in numbers:
                continue
            # Check equation 6: B + K = 130
            if B + K != 130:
                continue
            # Check equation 7: A > B
            if A <= B:
                continue
            # Check equation 11: A > K
            if A <= K:
                continue
            # If all conditions are satisfied, we have found a solution
            values = {'A': A, 'B': B, 'C': C, 'D': None, 'E': E, 'F': F, 'G': None, 'H': H, 'I': I, 'J': J, 'K': K}
            break
    if values:
        break

# Print the solution in alphabetical order
solution = [values[letter] for letter in sorted(values.keys())]
print(f"<<<{solution}>>>")