import math

def calculate_probability():
    """
    Calculates the probability that for each individual, there exists a type of 
    item for which they hold strictly more copies than any other individual.
    """

    # Use math.factorial for calculations
    fact = math.factorial

    # Step 1: Calculate S, the total number of possible distributions.
    # S = 25! / (5! * 5! * 5! * 5! * 5!)
    S = fact(25) // (fact(5)**5)

    # Step 2: Analyze and count favorable distributions for a fixed mapping.
    # We consider two main cases for the distribution matrix A.

    # Case 1: The matrix is 5*I (Identity).
    # A_ii = 5, A_ij = 0 for i!=j.
    # The number of ways to achieve this specific deal is:
    # ways(A) = (5!)^5 / (product of A_ij!). Here, (5!)^5 / (5!^5 * 0!^20) = 1.
    # There is only 1 such matrix for a fixed mapping.
    f_map_case1 = 1

    # Case 2: The matrix is 4*I + P, where P is a derangement permutation matrix.
    # The diagonal elements are 4, and each person has one item of another type.
    # The number of derangements of 5 items is D(5) = 44.
    # For each such matrix, the number of ways is:
    # ways(A) = (5!)^5 / (4!^5 * 1!^5) = (5)^5 = 3125.
    num_derangements_5 = 44
    ways_per_derangement_case = 5**5
    f_map_case2 = num_derangements_5 * ways_per_derangement_case

    # The sum of ways for a single, fixed specialization mapping (e.g., P_i specializes in T_i)
    f_map = f_map_case1 + f_map_case2
    
    # Step 3: Calculate F, the total number of favorable distributions.
    # We multiply by 5! for the number of possible specialization mappings.
    F = fact(5) * f_map
    
    # Step 4: Calculate the final probability P = F / S.
    # We can simplify the fraction by dividing by the greatest common divisor.
    common_divisor = math.gcd(F, S)
    F_simplified = F // common_divisor
    S_simplified = S // common_divisor
    
    # Final Probability
    P = F / S

    # Print the results as requested.
    print("Calculation Steps and Results:")
    print(f"- Total number of ways to distribute the items (S):")
    print(f"S = {S}")
    print("\n- Number of favorable distributions (F):")
    print(f"F = {F}")
    print("\n- The probability P is the ratio F / S:")
    print(f"P = {F} / {S}")
    print(f"   (Simplified: {F_simplified} / {S_simplified})")
    print(f"\nP ≈ {P}")

# Run the calculation and print the output.
calculate_probability()