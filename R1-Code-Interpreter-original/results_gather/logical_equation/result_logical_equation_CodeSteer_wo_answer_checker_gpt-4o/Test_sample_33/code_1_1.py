import itertools

def solve_constraints():
    numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]
    
    for perm in itertools.permutations(numbers):
        A, B, C, D, E, F, G, H, I, J, K = perm
        
        if (E - F == -59 and
            G - C == 6 and
            H + J == 107 and
            E + I == 137 and
            B == 4 * G and
            A + G == 234 and
            G == 3 * C and
            E - H == -64 and
            C - E == -13 and
            F + K == 123):
            return f"<<<{[A, B, C, D, E, F, G, H, I, J, K]}>>>"

print(solve_constraints())