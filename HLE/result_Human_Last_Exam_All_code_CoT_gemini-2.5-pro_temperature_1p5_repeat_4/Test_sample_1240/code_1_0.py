def solve_lattice_questions():
    """
    Analyzes the three questions about root systems of d-neighbors of Z^n.
    """

    # --- Question 1 ---
    # Is it true that for a d-neighbor N of Z^12, R_2(M) can be of type A_11?
    print("--- Analysis for Question 1 ---")
    n1 = 12
    root_system_type = "A_11"
    k1 = 11
    
    print(f"For n = {n1}, we ask if the visible root system R_2(M) can be {root_system_type}.")
    print("This requires the existence of a suitable integer d (the index of M in Z^n).")
    
    # A known theoretical condition is that if R_2(M) has an A_k component, then (k+1) must divide d.
    # Here k = 11.
    print(f"A necessary condition for an A_k component is that k+1 must divide d.")
    print(f"For A_{k1}, this means {k1}+1 = {k1+1} must divide d.")
    
    # We can choose d to satisfy this condition. Let's try d = 12.
    d1 = 12
    print(f"Let's test if a construction works for d = {d1}.")
    
    # Consider the lattice M = {v in Z^12 | sum of coordinates of v is 0 mod 12}.
    # The roots in Z^12 are of the form e_i - e_j or e_i + e_j.
    print("We construct M = {v in Z^12 | sum(v_i) = 0 mod 12}.")

    # Check which roots are in M.
    # For v = e_i - e_j, the sum of coordinates is 1 - 1 = 0.
    sum_em_ej = 0
    print(f"For a root v = e_i - e_j, the sum of coordinates is {sum_em_ej}.")
    print(f"Check: {sum_em_ej} mod {d1} = {sum_em_ej % d1}. This is 0, so these {k1*(k1+1)} roots (type A_11) are in M.")

    # For v = e_i + e_j, the sum of coordinates is 1 + 1 = 2.
    sum_ep_ej = 2
    print(f"For a root v = e_i + e_j, the sum of coordinates is {sum_ep_ej}.")
    print(f"Check: {sum_ep_ej} mod {d1} = {sum_ep_ej % d1}. This is not 0, so these roots are not in M.")
    
    print("Thus, R_2(M) is exactly of type A_11. This construction is valid.")
    answer_a = "Yes"
    print(f"Conclusion for (a): {answer_a}\n")

    # --- Question 2 ---
    # Can the visible root system R_2(M) of a d-neighbor N of Z^15 contain a D_7 component?
    print("--- Analysis for Question 2 ---")
    n2 = 15
    root_system_comp = "D_7"
    k2 = 7

    print(f"For n = {n2}, we ask if R_2(M) can contain a {root_system_comp} component for some d.")
    print("We can try to construct such a lattice M.")
    # Consider a lattice M built as a direct sum of the D_7 lattice and Z^8.
    # M = {(v_1,...,v_7, v_8,...,v_15) in Z^15 | sum of first 7 coordinates is even}
    # This corresponds to the lattice of the D_7 root system on the first 7 coordinates,
    # and the standard integer lattice Z^8 on the remaining 8.
    print(f"Consider the lattice M = (D_7 lattice)_coords(1-7) \u2295 (Z^8)_coords(8-15).")
    # The index d of this M in Z^15 is 2.
    d2 = 2
    print(f"This lattice M has index d = {d2} in Z^{n2}.")
    
    # Find the roots of this M. A vector v = (x, y) has norm^2 = x.x + y.y = 2.
    # Case 1: x.x=2, y.y=0. This means x is a root of D_7 and y is the zero vector. These are the D_7 roots.
    # Case 2: x.x=0, y.y=2. This means x is zero and y is a root of Z^8. The roots of Z^8 are of type D_8.
    print(f"The root system R_2(M) contains all vectors of the form (r, 0) where r is a D_7 root, and (0, s) where s is a D_8 root.")
    print("So, R_2(M) = D_7 \u2295 D_8.")
    print("This root system contains a D_7 component.")
    
    # This construction is known to correspond to a valid 2-neighbor.
    answer_b = "Yes"
    print(f"Conclusion for (b): {answer_b}\n")
    
    # --- Question 3 ---
    # For n=18 and d=5, is it possible for R_2(M) to include more than one D_n component?
    print("--- Analysis for Question 3 ---")
    n3 = 18
    d3 = 5
    
    print(f"For n = {n3} and d = {d3}, we ask if R_2(M) can have more than one D_k component.")
    # M is defined by a primitive glue vector u: M = {v in Z^18 | u.v = 0 mod 5}.
    print(f"The lattice M can be defined as M = {{v \u2208 Z^{n3} | u \u00B7 v \u2261 0 (mod {d3})}} for a primitive vector u.")
    
    # For a D_k component to exist on a set of coordinates I, its roots {e_a +- e_b | a,b in I} must be in M.
    print("For a D_k component on coordinates I, we need u_a \u00B1 u_b \u2261 0 (mod 5) for all a, b \u2208 I.")
    
    # This implies (u_a + u_b) + (u_a - u_b) = 2*u_a must be 0 mod 5.
    print("Adding (u_a + u_b) and (u_a - u_b), we get the equation: 2 \u00B7 u_a \u2261 0 (mod 5).")
    print("Since gcd(2, 5) = 1, this implies u_a \u2261 0 (mod 5).")
    
    print("So, for a D_k component on coordinates I, the glue vector u must have u_i \u2261 0 (mod 5) for all i \u2208 I.")
    
    # Now, assume there are two D_k components, on disjoint coordinate sets I_1 and I_2.
    print("If R_2(M) had two D_k components, on sets I_1 and I_2, then u_i \u2261 0 (mod 5) for all i \u2208 I_1 \u222A I_2.")
    
    # But if u_i = 0 mod 5 for all i in I_1 U I_2, then consider a root e_a + e_b where a is in I_1 and b is in I_2.
    # u . (e_a + e_b) = u_a + u_b = 0 + 0 = 0 mod 5.
    print("This implies that for any a \u2208 I_1 and b \u2208 I_2, the root e_a \u00B1 e_b is also in M.")
    print("This 'cross-root' connects the two supposed components, meaning they are not separate components but part of one larger D_k component.")
    
    print("Therefore, the root system R_2(M) can have at most one D_k component.")
    answer_c = "No"
    print(f"Conclusion for (c): {answer_c}\n")

    print("--- Final Answer ---")
    final_answer = f"(a) [{answer_a}]; (b) [{answer_b.lower()}]; (c) [{answer_c.lower()}]."
    print(final_answer)

solve_lattice_questions()
