def solve_csp():
    # Possible values for the letters
    values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

    # Iterate over possible values for C
    for C in values:
        if C % 2 != 0:
            continue  # C must be even

        # Calculate dependent variables
        D = C / 2.0
        if D not in values:
            continue

        K = 1.5 * C
        if K not in values:
            continue

        A = 3.6 * C
        if A not in values:
            continue

        L = 1.6 * C
        if L not in values:
            continue

        M = 4.8 * D
        if M not in values:
            continue

        F = 106 - C
        if F not in values:
            continue

        # Iterate over possible values for G
        for G in values:
            H = G + 42
            if H not in values:
                continue

            # Iterate over possible values for E
            for E in values:
                I = 57 - E
                if I not in values:
                    continue

                # Check if all values are unique
                if len({A, B, C, D, E, F, G, H, I, J, K, L, M}) == 13:
                    # All constraints are satisfied, print the result
                    print(f"<<<[{A}, {B}, {C}, {D}, {E}, {F}, {G}, {H}, {I}, {J}, {K}, {L}, {M}]>>>")
                    return

solve_csp()