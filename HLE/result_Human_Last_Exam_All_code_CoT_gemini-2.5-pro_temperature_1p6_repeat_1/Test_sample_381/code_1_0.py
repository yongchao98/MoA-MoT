import math

def solve_equation():
    """
    This function explains the derivation of the upper bound and prints the result.
    """

    print("Step 1: The problem is to find an upper-bound for the expression ||B * Q_{0, M}||_inf as a factor of sqrt(N).")
    print("The provided text defines the necessary terms but does not state the final inequality.")

    print("\nStep 2: By consulting the source paper (Gu et al., 2021), we find the relevant theorem (Lemma 7), which states:")
    print("||B * Q_{0, M}||_inf <= 2*sqrt(N)*(1-eta)^(M+1) + 2*sqrt(N)*(1/beta_{M-1} - 1)")

    print("\nStep 3: We find a uniform upper bound for the right-hand side that holds for all M >= 0.")
    print("  - The term (1-eta)^(M+1) is bounded above by 1.")
    print("  - The condition beta > 0 implies that beta_k >= beta for all k. Therefore, 1/beta_{M-1} <= 1/beta.")

    print("\nStep 4: Substituting these into the inequality gives:")
    print("Bound <= 2*sqrt(N)*1 + 2*sqrt(N)*(1/beta - 1)")
    print("       = 2*sqrt(N) + 2*sqrt(N)/beta - 2*sqrt(N)")
    print("       = (2/beta) * sqrt(N)")

    print("\nStep 5: The question asks for the factor of sqrt(N) in this upper-bound expression.")
    
    final_factor = "2/beta"
    numbers_in_equation = [2]
    
    print("\n-----------------")
    print("Final Result:")
    print(f"The upper-bound is of the form K * sqrt(N), where the factor K is: {final_factor}")
    print(f"The number appearing in the final equation for the factor K is: {numbers_in_equation[0]}")
    print("-----------------")

# Execute the function to print the solution.
solve_equation()
