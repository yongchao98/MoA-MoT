from itertools import permutations

numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

def check_conditions(A, B, C, D, E, F, G, H, I, J, K):
    return (H + J == 130 and
            C + I == 177 and  # Corrected this value
            C + G == 75 and
            abs(G - 3.0 * F) < 0.0001 and
            B - E == 150 and
            J - C == -18 and
            E + G == 123 and
            D + F == 96 and
            E + K == 78 and
            H + K == 124)

for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K = perm
    if check_conditions(A, B, C, D, E, F, G, H, I, J, K):
        print(f"<<<[{A}, {B}, {C}, {D}, {E}, {F}, {G}, {H}, {I}, {J}, {K}]>>>")
        break