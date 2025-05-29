# List of numbers to assign
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over possible values for J and K based on their relationship
for J in numbers:
    for K in numbers:
        if J == 1.4 * K:
            A = 4.0 * J
            if A in numbers:
                I = J - 4
                if I in numbers:
                    F = 57 - J
                    if F in numbers:
                        G = 60 - F
                        if G in numbers:
                            D = 34 - G
                            if D in numbers:
                                E = (J - 4) / 2.4
                                if E in numbers:
                                    H = 2.4 * E
                                    if H in numbers:
                                        C = 52 - H
                                        if C in numbers:
                                            L = 95 - F
                                            if L in numbers:
                                                B = next(x for x in numbers if x not in [A, B, C, D, E, F, G, H, I, J, K, L, M])
                                                M = next(x for x in numbers if x not in [A, B, C, D, E, F, G, H, I, J, K, L, B])
                                                print([A, B, C, D, E, F, G, H, I, J, K, L, M])
                                                break