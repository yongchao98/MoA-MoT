# Possible values for the letters
possible_values = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Iterate over possible values for C
for C in possible_values:
    J = 186 - C
    I = 157 - C
    B = C + 44
    K = C - 20
    D = 84 - C
    E = C - 27
    H = 261 - C
    G = C - 9
    F = (261 - C) / 3

    # Check if all values are in the possible values list and are unique
    values = [C, J, I, B, K, D, E, H, G, F]
    if all(v in possible_values for v in values) and len(set(values)) == len(values):
        # Check if F is an integer
        if F.is_integer():
            A = list(set(possible_values) - set(values))[0]  # The remaining value is A
            result = [A, B, C, D, E, F, G, H, I, J, K]
            print(result)
            break