def solve_lattice_questions():
    """
    This script analyzes three questions about root systems of d-neighbors of integer lattices.
    It explains the reasoning for each question and prints the final formatted answer.
    """

    print("--- Analyzing Question 1: Is R_2(M) of type A_11 possible for a d-neighbor of Z^12? ---")
    n1 = 12
    d1 = 3
    w1 = [1] * n1
    print(f"We test a construction for n={n1} with d={d1} and the glue vector w = {w1}.")
    print("The condition for a root v from Z^12 to be in the visible part R_2(M) is that w.v is divisible by d.")
    
    # Check roots of type e_i - e_j
    w_dot_v_minus = w1[0] - w1[1]
    print(f"For any root v = e_i - e_j, the dot product is w.v = w_i - w_j = {w1[0]} - {w1[1]} = {w_dot_v_minus}.")
    print(f"Checking divisibility by d={d1}: {w_dot_v_minus} % {d1} = {w_dot_v_minus % d1}. As this is 0, all A_11-type roots are included.")
    
    # Check roots of type e_i + e_j
    w_dot_v_plus = w1[0] + w1[1]
    print(f"For any root v = e_i + e_j, the dot product is w.v = w_i + w_j = {w1[0]} + {w1[1]} = {w_dot_v_plus}.")
    print(f"Checking divisibility by d={d1}: {w_dot_v_plus} % {d1} = {w_dot_v_plus % d1}. As this is not 0, these roots are excluded.")
    print("This construction results in a root system of type A_11. So the answer is Yes.")
    ans1 = "Yes"

    print("\n--- Analyzing Question 2: Can R_2(M) of a d-neighbor of Z^15 contain a D_7 component? ---")
    n2 = 15
    d2 = 2
    w2 = [1]*7 + [0]*8
    print(f"We test a construction for n={n2} with d={d2} and the glue vector w = {w2}.")
    print("To contain a D_7 component on the first 7 coordinates, all roots v = +/-e_i +/-e_j for 1<=i<j<=7 must satisfy w.v divisible by 2.")
    
    # Check D_7 roots on the first 7 coordinates, where w_i = 1
    w_i = 1
    w_j = 1
    w_dot_v_minus_d7 = w_i - w_j
    w_dot_v_plus_d7 = w_i + w_j
    
    print(f"For v = e_i - e_j (1<=i<j<=7), w.v = w_i - w_j = {w_i} - {w_j} = {w_dot_v_minus_d7}.")
    print(f"Checking divisibility: {w_dot_v_minus_d7} % {d2} = {w_dot_v_minus_d7 % d2}. Condition met.")
        
    print(f"For v = e_i + e_j (1<=i<j<=7), w.v = w_i + w_j = {w_i} + {w_j} = {w_dot_v_plus_d7}.")
    print(f"Checking divisibility: {w_dot_v_plus_d7} % {d2} = {w_dot_v_plus_d7 % d2}. Condition met.")
    print("Both conditions are met for the first 7 coordinates, so R_2(M) contains a D_7 component. The answer is yes.")
    ans2 = "yes"
    
    print("\n--- Analyzing Question 3: For n=18, d=5, can R_2(M) include more than one D_k component? ---")
    d3 = 5
    print(f"For n=18 and d={d3}, a D_k component on indices I requires w_i = c (mod 5) for all i in I, and 2*c must be divisible by 5.")
    
    c_sol = -1
    # Solve 2*c = 0 (mod 5)
    for c_val in range(d3):
        if (2 * c_val) % d3 == 0:
            c_sol = c_val
            break
            
    print(f"We solve the equation 2*c = 0 (mod 5) for c. The only solution is c = {c_sol}.")
    print("This means any D_k component must be on indices i where w_i is congruent to 0 (mod 5).")
    
    print("Now, assume we have two D_k components on disjoint index sets I1 and I2.")
    print("This requires w_i = 0 (mod 5) for i in I1, and w_j = 0 (mod 5) for j in I2.")
    print("For these components to be separate, the 'cross-term' roots (e.g., v = e_i +/- e_j for i in I1, j in I2) must be excluded.")
    
    w_i_comp = 0
    w_j_comp = 0
    cross_term_dot_sum = w_i_comp + w_j_comp
    print(f"Let's check the cross-term dot product: w.v = w_i + w_j = {w_i_comp} + {w_j_comp} = {cross_term_dot_sum} (mod 5).")
    
    if cross_term_dot_sum % d3 == 0:
        print(f"Since {cross_term_dot_sum} is divisible by {d3}, the cross-term roots are included, not excluded.")
        print("This connects the two components into a single larger D_{k1+k2} component.")
    
    print("Therefore, it is impossible to have more than one separate D_k component. The answer is no.")
    ans3 = "no"

    print("\n--------------------------------------------------")
    print("Final Answer Summary:")
    print(f"(a) [{ans1}]; (b) [{ans2}]; (c) [{ans3}].")
    print("<<< (a) [Yes]; (b) [yes]; (c) [no] >>>")

# Execute the function to solve the questions
solve_lattice_questions()