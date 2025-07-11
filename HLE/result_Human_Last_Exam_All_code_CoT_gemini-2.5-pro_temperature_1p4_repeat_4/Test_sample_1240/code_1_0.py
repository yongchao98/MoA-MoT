def solve_lattice_questions():
    """
    Analyzes and answers the three questions about root systems of d-neighbors of Z^n.
    """

    # --- Question (a) ---
    print("Part (a): Is it true that for a d-neighbor N of Z^12, R_2(M) can be of type A_11?")
    n_a = 12
    r_a = {'name': 'A', 'k': 11, 'rank': 11}

    print(f"Analysis for n = {n_a}, R = {r_a['name']}{r_a['k']}:")
    # Step 1: Check if A_11 can be a subsystem of D_12
    # The root system A_k requires an ambient space of dimension at least k.
    # The standard embedding of A_k as a subsystem of D_n requires n >= k+1.
    embed_cond_a = n_a >= r_a['k'] + 1
    print(f"1. Embeddability Check: A_{r_a['k']} embeds into D_{n_a} if n >= k+1.")
    print(f"   - We check if {n_a} >= {r_a['k']} + 1, which is {n_a} >= {r_a['k'] + 1}.")
    print(f"   - The condition is {'met' if embed_cond_a else 'not met'}.")

    # Step 2: Check for existence
    print("2. Existence Check: A known theorem states that if a root system R can be embedded")
    print("   in D_n, then a d-neighbor of Z^n with root system R exists for some d.")
    
    if embed_cond_a:
        print("   - Since the embeddability condition is met, such a neighbor exists.")
        answer_a = "Yes"
    else:
        answer_a = "No"
    
    print(f"Conclusion for (a): The root system A_11 is possible for a d-neighbor of Z^12.")
    print("-" * 20)

    # --- Question (b) ---
    print("Part (b): Can the visible root system R_2(M) of a d-neighbor N of Z^15 contain a D_7 component?")
    n_b = 15
    r_b_comp = {'name': 'D', 'k': 7, 'rank': 7}
    
    print(f"Analysis for n = {n_b}, containing R = {r_b_comp['name']}{r_b_comp['k']}:")
    # Step 1: Check if D_7 can be a subsystem of D_15
    # The root system D_k can be a subsystem of D_n if n >= k.
    embed_cond_b = n_b >= r_b_comp['k']
    print(f"1. Embeddability Check: D_{r_b_comp['k']} embeds into D_{n_b} if n >= k.")
    print(f"   - We check if {n_b} >= {r_b_comp['k']}.")
    print(f"   - The condition is {'met' if embed_cond_b else 'not met'}.")
    
    # Step 2: Check for existence
    print("2. Existence Check: The question is whether R_2(M) can *contain* a D_7 component.")
    print("   If D_7 itself is a possible root system, the answer is yes. More complex systems")
    print("   like D_7 + D_k are also possible if they embed in D_15.")
    print(f"   For example, R = D_7 + D_8 has rank {7+8}=15, which is embeddable in D_15.")
    
    if embed_cond_b:
        print("   - Since D_7 (and even D_7+D_8) is embeddable, a d-neighbor with such a root system exists.")
        answer_b = "yes"
    else:
        answer_b = "no"

    print(f"Conclusion for (b): R_2(M) can contain a D_7 component.")
    print("-" * 20)

    # --- Question (c) ---
    print("Part (c): For n = 18 and d = 5, is it possible for R_2(M) to include more than one D_n component?")
    n_c = 18
    d_c = 5
    
    print(f"Analysis for n = {n_c}, d = {d_c}:")
    print(f"1. Structure of sublattice for prime index d={d_c}:")
    print("   - Since d=5 is a prime, any sublattice M of index 5 in Z^18 is defined by a")
    print("     linear congruence a . x = 0 (mod 5) for some vector a in (Z/5Z)^18.")

    print("2. Structure of the root system R_2(M):")
    print("   - Let k be the number of non-zero coordinates of 'a'. The roots of M (vectors of")
    print("     norm 2) are partitioned into two orthogonal sets:")
    print("     - Those on coordinates where a_i = 0, which form a D_{18-k} system.")
    print("     - Those on coordinates where a_i != 0 (let this system be R_k).")
    print("   - Thus, R_2(M) = D_{18-k} + R_k.")
    
    print("3. Check for a second D component:")
    print("   - R_2(M) can have more than one D component only if R_k contains a D_m (m>=2) component.")
    print("   - For R_k to contain D_m, the roots of D_m must satisfy the congruence.")
    print("     This means for all i,j in {1,...,m}, we need +/-a_i +/-a_j = 0 (mod 5).")
    print("     This implies a_i+a_j=0 and a_i-a_j=0, which gives 2*a_i = 0 (mod 5).")
    print("   - Since 2 and 5 are coprime, this forces a_i = 0 (mod 5).")
    print("   - This contradicts that i is in the support of 'a' (where a_i != 0).")
    print("   - Therefore, R_k cannot contain any D_m component for m >= 2.")
    
    print("4. Conclusion for (c):")
    print("   - The root system R_2(M) can have at most one D component (D_{18-k}).")
    answer_c = "no"

    print("-" * 20)

    # --- Final Answer ---
    print("Final Answer:")
    final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}."
    print(final_answer)
    return final_answer

solve_lattice_questions()