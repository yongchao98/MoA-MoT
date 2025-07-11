import math

def solve_dice_problem():
    """
    Calculates the minimal expected number of rolls of a 7-sided die to
    generate a uniform random digit (0-9).
    """

    min_expected_rolls = float('inf')
    best_k = -1
    # Variables to store the parts of the fraction for the best k
    best_num = -1
    best_den = -1
    
    # Store all components for the final print statement
    equation_parts = {}

    print("Analyzing the expected number of rolls E(k) for different numbers of rolls per trial (k).")
    print("-" * 60)
    
    # We only need to check a few values of k, as E(k) tends to increase with k for k>=2.
    # k=1 is not possible as 7 < 10, so we start the search from k=2.
    for k in range(2, 6):
        base = 7
        target = 10
        
        # Total outcomes for k rolls
        total_outcomes = base**k
        
        # Number of outcomes we can map to the 10 digits uniformly
        accepted_outcomes = target * (total_outcomes // target)
        
        if accepted_outcomes == 0:
            print(f"For k={k}, total outcomes ({total_outcomes}) are less than target ({target}), so this is not a valid strategy.")
            continue

        # The expected number of rolls is k * (total_outcomes / accepted_outcomes)
        num = k * total_outcomes
        den = accepted_outcomes
        
        current_expected_rolls = num / den
        
        print(f"For k={k}:")
        print(f"  - Total outcomes = 7^{k} = {total_outcomes}")
        print(f"  - Accepted outcomes = 10 * floor({total_outcomes}/10) = {accepted_outcomes}")
        print(f"  - Expected rolls E({k}) = {k} * {total_outcomes} / {accepted_outcomes} = {num}/{den} ≈ {current_expected_rolls:.4f}")
        
        # Check if this k is the new minimum
        if current_expected_rolls < min_expected_rolls:
            min_expected_rolls = current_expected_rolls
            best_k = k
            best_num = num
            best_den = den
            equation_parts = {
                'k': k,
                'base': base,
                'target': target,
                'total_outcomes': total_outcomes,
                'accepted_outcomes': accepted_outcomes,
                'num': num,
                'den': den
            }

    # Simplify the fraction for the minimal expected value
    common_divisor = math.gcd(best_num, best_den)
    final_num = best_num // common_divisor
    final_den = best_den // common_divisor

    print("-" * 60)
    print(f"The minimal expected number of rolls is achieved with k = {best_k}.")
    print("\nThe calculation for the final answer is as follows:")
    
    # Unpack for clarity
    k_val = equation_parts['k']
    base_val = equation_parts['base']
    target_val = equation_parts['target']
    total_out = equation_parts['total_outcomes']
    accepted_out = equation_parts['accepted_outcomes']
    num_val = equation_parts['num']
    den_val = equation_parts['den']

    # Final formatted output of the equation as requested
    print(f"E = {k_val} * ({base_val}^{k_val}) / ({target_val} * floor({base_val}^{k_val} / {target_val}))")
    print(f"E = {k_val} * {total_out} / ({target_val} * {total_out // target_val})")
    print(f"E = {num_val} / {den_val}")
    print(f"E = {final_num}/{final_den}")

solve_dice_problem()