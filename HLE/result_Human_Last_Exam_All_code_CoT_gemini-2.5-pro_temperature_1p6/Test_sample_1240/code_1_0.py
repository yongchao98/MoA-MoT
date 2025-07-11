import numpy as np

def solve_questions():
    """
    Solves the three questions about root systems of d-neighbors of Z^n.
    The code provides step-by-step reasoning for each answer.
    """

    # A helper function to make print statements clearer
    def print_step(text):
        print(f"- {text}")

    def solve_q1():
        print("--- Question 1 ---")
        print("Is it true that for a d-neighbor N of Z^12, R_2(M) can be of type A_11?\n")
        
        n = 12
        d = 3
        u = np.ones(n, dtype=int)
        
        print_step("We test a potential construction.")
        print_step(f"Let n = {n}, d = {d}, and the primitive vector u = {list(u)}.")
        print_step(f"This defines the sublattice M = {{x in Z^{n} | u.x = 0 (mod {d})}}.")
        
        print("\nFirst, we check if a d-neighbor N can exist for this M.")
        u_dot_u = np.dot(u, u)
        # For even n, the condition is u.u = 0 mod (2d/gcd(2,d))
        mod = 2 * d // np.gcd(2, d)
        print_step(f"For n={n} (even), the existence condition is u.u = 0 (mod {mod}).")
        print_step(f"We calculate u.u = sum([1**2 for _ in range(12)]) = {u_dot_u}.")
        print_step(f"We check the condition: {u_dot_u} % {mod} = {u_dot_u % mod}.")
        if u_dot_u % mod == 0:
            print_step("The condition is satisfied, so such a neighbor N exists.")
        else:
            print_step("The condition is not satisfied. This construction is invalid.")
            print("\n(a) The answer based on this construction would be No.")
            return "No"
            
        print("\nNext, we determine the root system R_2(M).")
        print_step("The roots v in R_2(M) are vectors in Z^12 of squared norm 2 (i.e., v = +/- e_i +/- e_j) that also satisfy u.v = 0 (mod d).")

        # Check for A_11 type roots: e_i - e_j
        print_step(f"For roots of type v = e_i - e_j, u.v = u_i - u_j = {u[0]} - {u[1]} = {u[0] - u[1]}.")
        print_step(f"Since ({u[0]} - {u[1]}) % {d} = 0, all {n*(n-1)} vectors of the form e_i - e_j are in M. This forms an A_{n-1} = A_11 root system.")

        # Check for other roots like e_i + e_j
        print_step(f"For roots of type v = e_i + e_j, u.v = u_i + u_j = {u[0]} + {u[1]} = {u[0] + u[1]}.")
        print_step(f"Since ({u[0]} + {u[1]}) % {d} != 0, vectors of the form e_i + e_j are not in M.")

        print_step("Thus, the root system R_2(M) consists solely of the roots {+/- (e_i - e_j)}, which is precisely of type A_11.")
        print("\n(a) The answer is Yes.")
        return "Yes"

    def solve_q2():
        print("\n--- Question 2 ---")
        print("Can the visible root system R_2(M) of a d-neighbor N of Z^15 contain a D_7 component?\n")
        
        n = 15
        d = 4
        k = 7
        u = np.array([2]*k + [1]*(n-k))

        print_step("We test a potential construction.")
        print_step(f"Let n = {n}, d = {d}.")
        print_step(f"We partition the 15 coordinates into a set of {k} (S1) and a set of {n-k} (S2).")
        print_step(f"Let u_i = 2 for i in S1, and u_i = 1 for i in S2. So, u = {list(u)}.")
        is_primitive = np.gcd.reduce(u) == 1
        print_step(f"This u is primitive because gcd(2, 1) = {is_primitive}).")
        
        print("\nFirst, we check if a d-neighbor N can exist for this M.")
        u_dot_u = np.dot(u, u)
        # For odd n, the condition is u.u = 0 mod d
        mod = d
        print_step(f"For n={n} (odd), the existence condition is u.u = 0 (mod {mod}).")
        print_step(f"We calculate u.u = {k} * (2^2) + {n-k} * (1^2) = {k*4 + (n-k)} = {u_dot_u}.")
        print_step(f"We check the condition: {u_dot_u} % {mod} = {u_dot_u % mod}.")
        if u_dot_u % mod == 0:
            print_step("The condition is satisfied, so such a neighbor N exists.")
        else:
            print_step("The condition is not satisfied. This construction is invalid.")
            print("\n(b) The answer based on this construction would be No.")
            return "No"

        print("\nNext, we check for a D_7 component on S1.")
        print_step(f"A D_7 system requires vectors e_i +/- e_j to be in M for all i, j in S1.")
        i, j = 0, 1 # Example indices from S1
        
        uv_dot_minus = u[i] - u[j]
        print_step(f"For v = e_i - e_j (i,j in S1), u.v = u_i - u_j = {u[i]} - {u[j]} = {uv_dot_minus}.")
        print_step(f"Check: {uv_dot_minus} % {d} = {uv_dot_minus % d}. Condition holds.")
        
        uv_dot_plus = u[i] + u[j]
        print_step(f"For v = e_i + e_j (i,j in S1), u.v = u_i + u_j = {u[i]} + {u[j]} = {uv_dot_plus}.")
        print_step(f"Check: {uv_dot_plus} % {d} = {uv_dot_plus % d}. Condition holds.")
        print_step("So all vectors for a D_7 system are present.")

        print("\nFinally, we check if this D_7 system is a separate component.")
        print_step("This requires that no roots connect S1 and S2.")
        i_s1, j_s2 = 0, 7 # Example indices from S1 and S2
        
        uv_dot_cross_plus = u[i_s1] + u[j_s2]
        uv_dot_cross_minus = u[i_s1] - u[j_s2]
        print_step(f"For v = e_i +/- e_j (i in S1, j in S2), u.v = u_i +/- u_j = {u[i_s1]} +/- {u[j_s2]}.")
        print_step(f"This gives values {uv_dot_cross_plus} or {uv_dot_cross_minus}.")
        print_step(f"Neither {uv_dot_cross_plus} % {d} nor {uv_dot_cross_minus} % {d} is 0. So no connecting roots exist.")
        
        print_step("This construction yields a root system containing D_7 as a component.")
        print("\n(b) The answer is yes.")
        return "yes"
        
    def solve_q3():
        print("\n--- Question 3 ---")
        print("For n=18 and d=5, is it possible for R_2(M) to include more than one D_n component?\n")
        
        d = 5
        print_step("Let's analyze the necessary conditions for a D-type component.")
        print_step(f"Assume a D_k component exists on a set of coordinates S. For any i, j in S, both e_i - e_j and e_i + e_j must be in R_2(M).")
        print_step(f"This means u.(e_i - e_j) = 0 (mod {d}) and u.(e_i + e_j) = 0 (mod {d}).")
        print_step(f"This is equivalent to two conditions: (1) u_i - u_j = 0 (mod {d}) and (2) u_i + u_j = 0 (mod {d}).")
        print_step(f"Adding the two congruences gives 2*u_i = 0 (mod {d}).")
        print_step(f"Since gcd(2, {d}) = 1, we can divide by 2, which implies u_i = 0 (mod {d}).")
        print_step("Since this must hold for any i in S, it means all u_i for i in S must be multiples of 5.")

        print("\nNow, suppose we have two D-type components: D_k1 on S1 and D_k2 on S2.")
        print_step("From the above reasoning, we must have u_i = 0 (mod 5) for all i in S1, and u_j = 0 (mod 5) for all j in S2.")
        
        print("\nFor these two to be separate components, there should be no roots connecting S1 and S2.")
        print_step("A connecting root would be v = e_i +/- e_j, where i is in S1 and j is in S2.")
        print_step(f"Let's check the condition u.v = 0 (mod {d}) for such a root.")
        print_step(f"u.v = u_i +/- u_j. Since u_i = 0 (mod 5) and u_j = 0 (mod 5), their sum and difference are also 0 (mod 5).")
        print_step("For example, if u_i = 5a and u_j = 5b for integers a, b. Then u_i + u_j = 5*(a+b), so (u_i + u_j) % 5 = 0.")

        print_step("This means that connecting roots MUST exist. The supposed two components are actually connected into a single, larger D-type component.")
        print_step("Therefore, it is impossible for R_2(M) to have more than one D-type component for d=5.")

        print("\n(c) The answer is no.")
        return "no"

    a1 = solve_q1()
    a2 = solve_q2()
    a3 = solve_q3()

    print(f"\n<<<(a) {a1}; (b) {a2}; (c) {a3}>>>")

solve_questions()