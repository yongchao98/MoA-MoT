import math

def solve_diagonalisable_probability():
    """
    Calculates the probability that the given matrix is diagonalisable by numerically
    evaluating the infinite series derived from the problem's conditions.
    """
    # The terms in the series decrease very rapidly due to the factorial.
    # Summing up to 30 terms is more than sufficient for standard float precision.
    num_terms = 30

    # Calculate S1 = sum_{k=1 to inf} 1/(k!)^3
    sum_cub_gt0 = 0.0
    for k in range(1, num_terms):
        try:
            # Calculate term and add to sum
            term = 1.0 / (math.factorial(k)**3)
            sum_cub_gt0 += term
        except OverflowError:
            # Factorial grows very fast; stop if it exceeds float limits.
            break

    # Calculate S2 = sum_{k=0 to inf} 1/(k!)^2
    sum_sq = 0.0
    for k in range(num_terms):
        try:
            # Calculate term and add to sum
            term = 1.0 / (math.factorial(k)**2)
            sum_sq += term
        except OverflowError:
            break
            
    e_minus_1 = math.exp(-1)

    # Calculate P(not diagonalisable) using the derived formula
    prob_not_diag = e_minus_1 * sum_cub_gt0 / sum_sq

    # The probability of being diagonalisable is 1 minus the probability of not being diagonalisable
    prob_diag = 1 - prob_not_diag

    print("The probability is calculated by the formula:")
    print("P(diagonalisable) = 1 - e^(-1) * (sum_{k=1 to inf} 1/(k!)^3) / (sum_{k=0 to inf} 1/(k!)^2)")
    print("\nSubstituting the computed values for the components:")
    print(f"e^(-1) = {e_minus_1}")
    print(f"sum_{k=1 to inf} 1/(k!)^3 ≈ {sum_cub_gt0}")
    print(f"sum_{k=0 to inf} 1/(k!)^2 ≈ {sum_sq}")
    
    print("\nThe final equation with the computed numbers is:")
    print(f"P(diagonalisable) = 1 - {e_minus_1} * {sum_cub_gt0} / {sum_sq}")
    
    print("\nResult:")
    print(f"The probability that the matrix is diagonalisable is: {prob_diag}")

solve_diagonalisable_probability()