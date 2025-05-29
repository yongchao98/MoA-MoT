from sympy import symbols, Eq, solve
from itertools import permutations

# Define symbols
H, I, E, K = symbols('H I E K')

# Solve H + I = 46 and I = 3.6H
eq1 = Eq(H + I, 46)
eq2 = Eq(I, 3.6 * H)
solutions_HI = solve((eq1, eq2), (H, I))

# Solve K = 3.2E
eq3 = Eq(K, 3.2 * E)
solutions_KE = solve(eq3, K)

# Available numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Filter possible values for H, I
possible_HI = []
for sol in solutions_HI:
    if isinstance(sol, tuple):
        h, i = sol
        if h.is_real and i.is_real and h in numbers and i in numbers:
            possible_HI.append((int(h), int(i)))

# Filter possible values for E, K
possible_KE = []
for e in numbers:
    k = solve(eq3.subs(E, e), K)
    if k and isinstance(k, list) and k[0].is_real and k[0] in numbers:
        possible_KE.append((e, int(k[0])))

# Iterate over possible values
for h, i in possible_HI:
    for e, k in possible_KE:
        if k == 3.2 * e:
            # Remaining numbers after assigning H, I, E, K
            remaining_numbers = [n for n in numbers if n not in {h, i, e, k}]
            
            # Iterate over permutations of remaining numbers
            for perm in permutations(remaining_numbers):
                A, B, C, D, F, G, J, L, M = perm
                
                # Check all conditions
                if (B + L == 35 and
                    J > D and
                    L - J == -8 and
                    B + I == 64 and
                    C + L == 31 and
                    G + K == 19 and
                    M - J == 30):
                    
                    # If all conditions are satisfied, print the result
                    print(f"<<<[{A}, {B}, {C}, {D}, {e}, {F}, {G}, {h}, {i}, {J}, {k}, {L}, {M}]>>>")
                    break