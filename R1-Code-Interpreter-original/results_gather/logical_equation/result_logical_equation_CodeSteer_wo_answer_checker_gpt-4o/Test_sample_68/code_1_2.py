from sympy import symbols, Eq, solve

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
possible_HI = [(h.evalf(), i.evalf()) for h, i in solutions_HI if h.evalf() in numbers and i.evalf() in numbers]

# Filter possible values for E, K
possible_KE = [(e, k.evalf()) for e in numbers for k in solutions_KE if e in numbers and k.evalf() in numbers]

# Iterate over possible values
for h, i in possible_HI:
    for e, k in possible_KE:
        if k == 3.2 * e:
            # Remaining numbers after assigning H, I, E, K
            remaining_numbers = [n for n in numbers if n not in {h, i, e, k}]
            
            # Iterate over permutations of remaining numbers
            from itertools import permutations
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