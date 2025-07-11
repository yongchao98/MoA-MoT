import math

def solve_probability():
    """
    Calculates the probability based on the simplifying assumption explained above.
    """
    
    # Step 1: Calculate S, the total number of ways to distribute the items.
    # S = 25! / (5!)^5
    try:
        S = math.factorial(25) // (math.factorial(5)**5)
    except OverflowError:
        print("Error: Calculation involves numbers too large to handle directly.")
        return

    # Step 2: Calculate F, the number of favorable distributions based on our assumption.
    
    # The number of ways to assign dominant types to individuals is 5!
    num_assignments = math.factorial(5)
    
    # Calculate F_fixed, the number of favorable ways for one fixed assignment.
    # It's the sum of ways for the two assumed cases.
    
    # Case 1: The "perfect" distribution (A=5I).
    # The number of ways to form this specific deal is 1.
    W_case1 = 1
    
    # Case 2: The "near-perfect" distribution (A = 4I + P_d).
    # The number of ways to form one such deal is 5^5.
    W_case2_per_matrix = 5**5
    
    # The number of derangements of 5 items.
    D5 = 44
    
    # Total ways for Case 2 is D5 * W_case2_per_matrix.
    Total_W_case2 = D5 * W_case2_per_matrix
    
    # F_fixed is the sum of ways for both cases.
    F_fixed = W_case1 + Total_W_case2
    
    # The total number of favorable distributions, F.
    F = num_assignments * F_fixed
    
    # Step 3: Calculate and print the final probability P = F / S.
    print("This solution is based on a strong simplifying assumption about the nature of favorable distributions.")
    print(f"Total number of possible distributions, S = {S}")
    print(f"Number of favorable distributions under the assumption, F = {F}")
    print("\nThe final probability P is the ratio F / S:")
    print(f"P = {F} / {S}")
    
    # The final probability value
    P = F / S
    print(f"P ≈ {P}")

solve_probability()