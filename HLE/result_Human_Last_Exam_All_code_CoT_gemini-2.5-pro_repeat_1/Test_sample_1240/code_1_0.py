def solve_lattice_questions():
    """
    Solves and explains the answers to the three questions about root systems
    of d-neighbors of Z^n.
    """

    # --- Question 1 ---
    print("--- Question 1 ---")
    print("Is it true that for a d-neighbor N of Z^12, R_2(M) can be of type A_11?")
    print("Plan: We will try to construct a suitable lattice M = {x in Z^12 | x.u = 0 mod d} whose root system is A_11.\n")

    n1 = 12
    d1 = 3
    u1 = [1] * n1
    
    print(f"Let n = {n1}. We propose a configuration with d = {d1} and the primitive vector u = {tuple(u1)}.")
    print("A root v = a*e_i + b*e_j is in R_2(M) if a*u_i + b*u_j is divisible by d.")
    print(f"With our choice, the condition is: a*{u1[0]} + b*{u1[1]} = a + b must be divisible by {d1}.\n")

    print("Analysis of root types:")
    # Case 1: A_11 roots (e.g., e_i - e_j)
    a, b = 1, -1
    result1 = a + b
    print(f"1. For A_11-type roots like e_i - e_j, we have (a,b) = (1,-1) or (-1,1).")
    print(f"   The condition check is: a + b = {a} + ({b}) = {result1}.")
    print(f"   Since {result1} is divisible by {d1}, these roots are INCLUDED.")

    # Case 2: Other roots (e.g., e_i + e_j)
    a, b = 1, 1
    result2 = a + b
    print(f"2. For other roots like e_i + e_j, we have (a,b) = (1,1) or (-1,-1).")
    print(f"   The condition check is: a + b = {a} + {b} = {result2}.")
    print(f"   Since {result2} is not divisible by {d1}, these roots are EXCLUDED.")
    
    print("\nConclusion: This construction yields exactly the root system A_11.")
    answer_a = "Yes"
    print(f"(a) The answer is {answer_a}.")

    # --- Question 2 ---
    print("\n\n--- Question 2 ---")
    print("Can the visible root system R_2(M) of a d-neighbor N of Z^15 contain a D_7 component?")
    print("Plan: We will construct a suitable M that has a D_7 component.\n")
    
    n2 = 15
    k2 = 7 # For D_7
    d2 = 3
    u2 = [0] * k2 + [1] * (n2 - k2)

    print(f"Let n = {n2}. We want a D_{k2} component. A D_k component on an index set I requires that 2*u_i is divisible by d for all i in I.")
    print(f"Let's choose d = {d2}. The condition 2*u_i = 0 mod {d2} implies u_i = 0 mod {d2}.")
    print(f"We define a primitive vector u by setting u_i=0 for the first {k2} indices and u_i=1 for the remaining {n2-k2}.")
    print(f"u = {tuple(u2)}\n")

    print("Analysis of root connections:")
    # Case 1: Roots within the first 7 indices
    ui, uj = 0, 0
    print(f"1. For roots with indices i,j <= {k2}, the condition is a*u_i + b*u_j = a*{ui} + b*{uj} = 0.")
    print(f"   This is always divisible by {d2}, so all vectors +/-e_i +/-e_j are roots. This forms a D_7 system.")

    # Case 2: Roots connecting the two sets of indices
    ui, uj = 0, 1
    result_conn = 1*ui + 1*uj
    print(f"2. For roots connecting i <= {k2} and j > {k2}, the condition is a*u_i + b*u_j = a*{ui} + b*{uj} = +/-{result_conn}.")
    print(f"   This is not divisible by {d2}, so there are no connecting roots. The D_7 system is a separate component.")

    print("\nConclusion: We have successfully constructed an R_2(M) containing a D_7 component.")
    answer_b = "yes"
    print(f"(b) The answer is {answer_b}.")

    # --- Question 3 ---
    print("\n\n--- Question 3 ---")
    print("For n = 18 and d = 5, is it possible for R_2(M) to include more than one D_n component?")
    print("We interpret this as asking if R_2(M) can have more than one D-type component (e.g., D_k and D_l).\n")
    
    n3 = 18
    d3 = 5
    print(f"Let n = {n3} and d = {d3}. A D_k component on an index set I requires 2*c = 0 mod d, where u_i = c mod d for all i in I.")
    print(f"The equation we must solve is 2 * c = 0 mod {d3}.")
    c_val = 0
    print(f"Since gcd(2, {d3}) = 1, the only solution is c = {c_val} mod {d3}.")
    print("This means that any D-type component must be on a set of indices I where u_i = 0 mod 5 for all i in I.\n")

    print("Proof by contradiction:")
    print("1. Assume there are two separate D-components, on disjoint index sets I and J.")
    print("2. This implies u_i = 0 mod 5 for all i in I, and u_j = 0 mod 5 for all j in J.")
    print("3. Let's check for a root v connecting I and J, for instance v = e_i + e_j with i in I, j in J.")
    ui, uj = 0, 0
    print(f"   The condition is v.u = u_i + u_j = {ui} + {uj} = 0.")
    print("4. Since 0 is divisible by 5, this connecting root must exist in R_2(M).")
    print("5. The existence of a connecting root means I and J are part of the same component, which contradicts our assumption.")
    
    print("\nConclusion: The set of all indices k where u_k = 0 mod 5 forms a single D-type component. It's impossible to have more than one.")
    answer_c = "no"
    print(f"(c) The answer is {answer_c}.")

    # --- Final Formatted Answer ---
    final_answer_str = f"(a) [{answer_a}]; (b) [{answer_b}]; (c) [{answer_c}]."
    print(f"\n<<<{final_answer_str}>>>")

solve_lattice_questions()