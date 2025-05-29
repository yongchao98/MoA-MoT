possible_values = [3, 6, 9, 20, 32, 36, 80, 120, 121]

for H in possible_values:
    I = 241 - H
    B = H - 88
    D = 141 - I
    A = 201 - I
    G = H - 111
    C = G / 1.5
    E = 1.8 * D
    
    # Check if all values are integers and in the possible_values list
    values = [A, B, C, D, E, G, H, I]
    if all(v in possible_values for v in values) and len(set(values)) == 8:
        print(f"<<<{[A, B, C, D, E, G, H, I]}>>>")
        break