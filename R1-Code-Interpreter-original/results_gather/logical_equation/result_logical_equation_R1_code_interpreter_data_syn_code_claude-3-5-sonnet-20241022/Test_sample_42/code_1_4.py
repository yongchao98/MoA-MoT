from itertools import permutations

numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K = perm
    if (H + J == 130 and
        C + I == 177 and
        C + G == 75 and
        abs(G - 3.0 * F) < 0.0001 and
        B - E == 150 and
        J - C == -18 and
        E + G == 123 and
        D + F == 96 and
        E + K == 78 and
        H + K == 124):
        print(f"<<<[{A}, {B}, {C}, {D}, {E}, {F}, {G}, {H}, {I}, {J}, {K}]>>>")
        break