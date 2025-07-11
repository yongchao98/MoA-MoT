def solve_case1():
    """
    Solves for the maximum N_r in Case 1, where phi(infinity) = 0.
    The constraints are m * N_r = 6 and |V_1| = N_r * (m - 1) - 1 >= 0.
    """
    print("--- Analyzing Case 1: phi(infinity) = 0 ---")
    print("The degree equation gives m * N_r = 6.")
    
    max_Nr = 0
    # m must be a divisor of 6
    divisors_of_6 = [1, 2, 3, 6]
    
    for m in divisors_of_6:
        Nr = 6 // m
        print(f"Checking m = {m}, which implies N_r = {Nr}.")
        
        # Check the constraint from the Riemann-Hurwitz formula
        num_q_vertices = Nr * (m - 1) - 1
        
        if num_q_vertices >= 0:
            print(f"  VALID: Number of q-vertices |V_1| = {num_q_vertices}, which is non-negative.")
            if Nr > max_Nr:
                max_Nr = Nr
        else:
            print(f"  INVALID: Number of q-vertices |V_1| = {num_q_vertices}, which is negative.")
            
    print(f"\nMaximum N_r found in Case 1 is: {max_Nr}")
    return max_Nr

def solve_case2():
    """
    Solves for the maximum N_r in Case 2, where phi(infinity) = 1.
    The constraints are m * N_r = 4 and |V_1| = N_r * (m - 1) >= 1.
    """
    print("\n--- Analyzing Case 2: phi(infinity) = 1 ---")
    print("The degree equation gives m * N_r = 4.")
    
    max_Nr = 0
    # m must be a divisor of 4
    divisors_of_4 = [1, 2, 4]
    
    for m in divisors_of_4:
        Nr = 4 // m
        print(f"Checking m = {m}, which implies N_r = {Nr}.")
        
        # Check the constraint from the Riemann-Hurwitz formula
        # We need |V_1| >= 1 because infinity itself is a q-vertex.
        num_q_vertices = Nr * (m - 1)
        
        if num_q_vertices >= 1:
            print(f"  VALID: Number of q-vertices |V_1| = {num_q_vertices}, which is >= 1.")
            if Nr > max_Nr:
                max_Nr = Nr
        else:
            print(f"  INVALID: Number of q-vertices |V_1| = {num_q_vertices}, which is less than 1.")
    
    print(f"\nMaximum N_r found in Case 2 is: {max_Nr}")
    return max_Nr

if __name__ == '__main__':
    max_Nr_case1 = solve_case1()
    max_Nr_case2 = solve_case2()
    
    overall_max = max(max_Nr_case1, max_Nr_case2)
    
    print("\n--- Conclusion ---")
    print(f"Comparing the maximum N_r from both cases ({max_Nr_case1} and {max_Nr_case2}).")
    print(f"The overall maximum number of vertices labelled 'r' within ]0, 1[ is {overall_max}.")
    
    # Final Answer
    # The final equation is m(N_r - N_p - n_p) = 6 or 4.
    # With m=2, Nr=3 we get 2*(3-0-0) = 6.
    print("\nThe final answer is derived from the equation m * N_r = 6, with m=2, leading to N_r=3.")
    print("Final Equation Check: m=2, N_r=3. Number of q-vertices |V_1| = 3*(2-1)-1=2, which is valid.")
    print("The final equation is 2 * 3 = 6")
