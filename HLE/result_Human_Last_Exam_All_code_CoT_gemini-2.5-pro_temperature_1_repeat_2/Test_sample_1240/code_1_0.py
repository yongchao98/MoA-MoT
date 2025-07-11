import math

def is_perfect_square(n):
    """Checks if a number is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = int(math.sqrt(n))
    return x * x == n

def solve_and_print():
    """Solves the three questions and prints the results with reasoning."""
    
    # --- Question 1 ---
    # Is it true that for a d-neighbor N of Z^12, R2(M) can be of type A_11?
    n_a = 11
    det_A_n = n_a + 1
    if is_perfect_square(det_A_n):
        ans1 = "Yes"
    else:
        ans1 = "No"

    # --- Question 2 ---
    # Can the visible root system R2(M) of a d-neighbor N of Z^15 contain a D_7 component?
    # Conditions: 1. Embedding possible (rank, det). 2. Glue theory holds for some d.
    # rank(D_7) = 7 <= 15. OK. det(D_7) = 4 (perfect square). OK.
    # Glue theory: exp(D_7*/D_7) | d. exp is 4. So we need 4|d.
    # Since a d=4 exists, it is possible.
    ans2 = "Yes"
    
    # --- Question 3 ---
    # For n = 18 and d = 5, is it possible for R2(M) to include more than one D_n component?
    # Condition: exp(D_k*/D_k) | d for any k. d=5.
    d_q3 = 5
    exp_D_even = 2
    exp_D_odd = 4
    if (d_q3 % exp_D_even == 0) or (d_q3 % exp_D_odd == 0):
        ans3 = "Yes" 
    else:
        # Neither 2 nor 4 divides 5, so no D_k component is possible at all.
        ans3 = "No"
        
    # Print the step-by-step reasoning
    print("(a) To determine if R2(M) can be of type A_11, we check if the A_11 root lattice can be a sublattice of Z^12.")
    print("A necessary condition for a lattice to be a sublattice of Z^n is that its determinant must be a perfect square.")
    print(f"The determinant of an A_n root lattice is calculated as n + 1.")
    print(f"For A_11, the determinant is {n_a} + 1 = {det_A_n}.")
    print(f"We check if {det_A_n} is a perfect square: {is_perfect_square(det_A_n)}.")
    print("Since the determinant is not a perfect square, the A_11 lattice cannot be a sublattice of Z^12. Therefore, R2(M) cannot be of type A_11.")
    
    print("\n(b) To determine if R2(M) can contain a D_7 component, we check two conditions.")
    print("First, the D_7 lattice must be embeddable in Z^15. Its rank is 7 (<= 15), and its determinant is 4, which is a perfect square. So, embedding is possible.")
    print("Second, according to 'glue theory', a d must exist such that d is divisible by the exponent of the discriminant group of the D_7 lattice.")
    print("For a D_k lattice with k odd (here k=7), the exponent is 4.")
    print(f"We need to find an integer d such that 4 divides d. A valid choice is d = 4, since 4 % 4 = 0.")
    print("Since a suitable d exists for which the necessary conditions are met, the answer is yes.")

    print("\n(c) To determine if R2(M) can include D_k components for d = 5, we use the same 'glue theory' condition.")
    print(f"The condition is that the exponent of the D_k lattice's discriminant group must divide d = {d_q3}.")
    print(f"For a D_k component with k even, the exponent is {exp_D_even}. Does {exp_D_even} divide {d_q3}? {d_q3 % exp_D_even == 0}.")
    print(f"For a D_k component with k odd, the exponent is {exp_D_odd}. Does {exp_D_odd} divide {d_q3}? {d_q3 % exp_D_odd == 0}.")
    print("Since neither 2 nor 4 divides 5, no D_k component can exist in the root system for d=5.")
    print("Therefore, it is impossible for R2(M) to include more than one (or even one) D_k component.")

    final_answer_str = f"(a) [{ans1}]; (b) [{ans2.lower()}]; (c) [{ans3.lower()}]."
    print("\n--------------------------------")
    print("Final Answer:")
    print(final_answer_str)
    print(f"<<<{final_answer_str}>>>")

solve_and_print()