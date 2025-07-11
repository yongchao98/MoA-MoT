import math

def solve_distribution_probability():
    """
    This function calculates the probability that for each of five individuals,
    there exists a unique type of item for which they hold strictly more copies
    than any other individual.

    The problem involves distributing 25 items (5 types, 5 copies each)
    among 5 individuals (5 items each).
    """

    # Step 1: Calculate S, the total number of possible distributions.
    # This is the number of ways to arrange 25 items with 5 groups of 5 identical items.
    # S = 25! / (5! * 5! * 5! * 5! * 5!) = 25! / (5!)^5
    try:
        s_val = math.factorial(25) // (math.factorial(5)**5)
    except OverflowError:
        print("Error: Calculation involves numbers too large to handle.")
        return

    # Step 2: Calculate F, the number of favorable distributions.
    # This requires counting distribution matrices M that satisfy the condition.
    # We sum the number of permutations for each valid matrix structure,
    # then multiply by 5! for the choice of specialty types.

    # Number of permutations for a given matrix M is (5!)^5 / product(M_ij!)
    
    # Contribution from Type 1 matrices (M = 5I): 1 matrix
    # N_perms = 1
    n_perms_sum = 1 * 1
    
    # Contribution from Type 2 matrices (M = 4I + P_derangement): !5 = 44 matrices
    # N_perms = 5^5
    n_perms_sum += 44 * (5**5)
    
    # Contribution from Type 3 matrices (M = 3I + 2*P_5-cycle): (5-1)! = 24 matrices
    # N_perms = 10^5
    n_perms_sum += 24 * (10**5)
    
    # Contribution from Type 4 matrices (M = 3I + Adj(C5)): (5-1)!/2 = 12 matrices
    # N_perms = 20^5
    n_perms_sum += 12 * (20**5)
    
    # Total F is 5! times this sum
    f_val = math.factorial(5) * n_perms_sum
    
    # Step 3: Calculate the probability P = F / S
    probability = f_val / s_val
    
    # Output the components of the final equation as requested.
    print(f"F = {f_val}")
    print(f"S = {s_val}")
    print(f"P = F / S = {f_val} / {s_val} = {probability}")

solve_distribution_probability()