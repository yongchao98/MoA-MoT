import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return abs(a)

def gcd_list(numbers):
    """Computes the greatest common divisor of a list of numbers."""
    if not numbers:
        return 0
    result = numbers[0]
    for i in range(1, len(numbers)):
        result = gcd(result, numbers[i])
    return result

def solve_lattice_questions():
    """
    Analyzes and answers the three questions about root systems of d-neighbors.
    """
    
    # --- Part (a) ---
    print("Part (a): Is it true that for a d-neighbor N of Z^12, R_2(M) can be of type A_11?")
    
    n_a = 12
    ans_a = "No"

    # We need to find a primitive vector u and integer d such that:
    # 1. u.u = 0 (mod d)
    # 2. For roots v = e_i - e_j (A_11), v.u = u_i - u_j = 0 (mod d)
    # 3. For roots v = e_i + e_j, v.u = u_i + u_j != 0 (mod d)
    #
    # Condition (2) implies u_i = u_j (mod d) for all i,j.
    # A simple choice is to have all components of u be equal. For primitivity, let's try u = (1, 1, ..., 1).
    u_a = [1] * n_a
    u_dot_u_a = sum(x * x for x in u_a)
    
    # Check primitivity of u
    if gcd_list(u_a) == 1:
        # Condition (1) requires u.u = 12 = 0 (mod d). Let's choose d = 12.
        d_a = 12
        if u_dot_u_a % d_a == 0:
            # Check root conditions
            dot_A11 = u_a[0] - u_a[1]  # Represents u_i - u_j
            dot_D12_extra = u_a[0] + u_a[1] # Represents u_i + u_j

            if dot_A11 % d_a == 0 and dot_D12_extra % d_a != 0:
                ans_a = "Yes"
                print(f"Let n = {n_a}. We choose the characteristic vector u = {tuple(u_a)} and d = {d_a}.")
                print(f"This vector u is primitive (gcd = {gcd_list(u_a)}).")
                print(f"The condition u.u = 0 (mod d) becomes {u_dot_u_a} = 0 (mod {d_a}), which is true since {u_dot_u_a} % {d_a} = {u_dot_u_a % d_a}.")
                print("The roots of A_11 are v = e_i - e_j. For these, v.u = u_i - u_j.")
                print(f"  v.u = {u_a[0]} - {u_a[1]} = {dot_A11}. Since {dot_A11} mod {d_a} is 0, all A_11 roots are included.")
                print("The other roots in D_12 are v = e_i + e_j. For these, v.u = u_i + u_j.")
                print(f"  v.u = {u_a[0]} + {u_a[1]} = {dot_D12_extra}. Since {dot_D12_extra} mod {d_a} is not 0, these roots are excluded.")
                print("Thus, with this construction, R_2(M) is precisely of type A_11.")
    print(f"Result for (a): {ans_a}\n")


    # --- Part (b) ---
    print("Part (b): Can the visible root system R_2(M) of a d-neighbor of Z^15 contain a D_7 component?")
    n_b = 15
    ans_b = "No"

    # A D_k component on indices I requires u_i +/- u_j = 0 (mod d) for all i,j in I.
    # This implies u_i = c (mod d) and 2c = 0 (mod d).
    # To satisfy 2c = 0 (mod d) for c not being 0 (mod d), d must be even. Let's try d = 2.
    d_b = 2
    # The condition 2c = 0 (mod 2) holds for any c. Let's take c=1.
    # So we need u_i = 1 (mod 2) for the 7 indices of the D_7 component.
    # We also need u.u = 0 (mod 2), meaning an even number of odd components in u.
    # With 7 odd components already required, we need at least one more. Let's have 8.
    num_ones = 8
    u_b = [1] * num_ones + [0] * (n_b - num_ones)
    u_dot_u_b = sum(x * x for x in u_b)

    if gcd_list(u_b) == 1 and u_dot_u_b % d_b == 0:
        # Check if D_7 is formed for indices 1..7
        # For i,j in {1..7}, u_i = 1, u_j = 1.
        dot_sum = 1 + 1
        dot_diff = 1 - 1
        if dot_sum % d_b == 0 and dot_diff % d_b == 0:
            ans_b = "yes"
            print(f"Let n = {n_b}. We need to find a suitable (u, d).")
            print(f"A D_k component requires 2c = 0 (mod d). Let's choose d = {d_b}, c = 1.")
            print(f"The condition u.u = 0 (mod {d_b}) requires an even number of odd components in u.")
            print("For a D_7 component, we need 7 indices i with u_i=1 (mod 2). We add one more to make the count of odd components 8 (even).")
            print(f"Let u = {tuple(u_b)}. This vector is primitive (gcd = {gcd_list(u_b)}).")
            print(f"u.u = {u_dot_u_b}, and {u_dot_u_b} = 0 (mod {d_b}) is true.")
            print("For any i, j in the first 7 indices, u_i = 1, u_j = 1. The dot products are:")
            print(f"  u_i + u_j = {1} + {1} = {dot_sum}, which is 0 (mod {d_b})")
            print(f"  u_i - u_j = {1} - {1} = {dot_diff}, which is 0 (mod {d_b})")
            print("Both conditions hold, so a D_7 component exists.")
    print(f"Result for (b): {ans_b}\n")
    

    # --- Part (c) ---
    print("Part (c): For n = 18 and d = 5, is it possible for R_2(M) to include more than one D_n component?")
    n_c = 18
    d_c = 5
    ans_c = "no"
    
    # A D_k component requires u_i=c (mod d) and 2c = 0 (mod d).
    print(f"Let n = {n_c} and d = {d_c}. A D_k component on indices I requires u_i = c (mod {d_c}) for i in I and 2c = 0 (mod {d_c}).")
    print(f"The equation 2c = 0 (mod 5) has a unique solution because gcd(2, 5) = 1.")
    print("Multiplying by 3 (the inverse of 2 mod 5), we get 6c = 0 (mod 5), which simplifies to c = 0 (mod 5).")
    print("This means that any set of indices supporting a D-type component must have u_i = 0 (mod 5).")
    print("Suppose we have two D-components on disjoint index sets P1 and P2.")
    print("For any index i in P1, u_i = 0 (mod 5). For any index j in P2, u_j = 0 (mod 5).")
    print("Now consider a 'cross-root' v = e_i + e_j. Its dot product with u is v.u = u_i + u_j.")
    print(f"v.u = {0} + {0} = {0}. This is 0 (mod {d_c}).")
    print("This means that roots exist connecting P1 and P2, so they are not separate components.")
    print("All indices with u_i = 0 (mod 5) form a single, large D-type component. Thus, more than one is impossible.")
    print(f"Result for (c): {ans_c}\n")

    final_answer = f"(a) [{ans_a}]; (b) [{ans_b}]; (c) [{ans_c}]."
    print(f"<<<{final_answer}>>>")

solve_lattice_questions()