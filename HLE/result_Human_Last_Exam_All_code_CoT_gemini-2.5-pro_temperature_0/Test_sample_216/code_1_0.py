import math

def calculate_performance_bound(H, lmbda):
    """
    Calculates the tightest upper bound on the performance difference J(pi*) - J(pi_hat).

    This calculation is based on the following derivation:
    1. The performance difference in imitation learning is bounded by the average
       policy error (epsilon_TV):
       J(pi*) - J(pi_hat) <= 2 * H^2 * epsilon_TV
       (assuming rewards are normalized to a maximum of 1).

    2. The problem states that the population total variation (TV) risk T
       is bounded: T <= |A| * (1 - exp(-lambda)).

    3. We interpret the non-standard risk T as T = |A| * epsilon_TV. This makes
       the problem statement consistent, as it implies:
       |A| * epsilon_TV <= |A| * (1 - exp(-lambda))
       which simplifies to:
       epsilon_TV <= 1 - exp(-lambda).

    4. Substituting this into the performance bound gives the final result:
       J(pi*) - J(pi_hat) <= 2 * H^2 * (1 - exp(-lambda)).

    Args:
        H (int): The episode horizon.
        lmbda (float): The hyperparameter lambda from the algorithm.
    """
    if H < 0 or lmbda < 0:
        print("Error: H and lambda must be non-negative.")
        return

    # The final formula for the upper bound
    # Note: The size of the action space |A| cancels out in the derivation.
    bound = 2 * (H**2) * (1 - math.exp(-lmbda))

    # Output the formula and the calculation as requested
    print("The formula for the tightest upper bound is:")
    print("J(pi*) - J(pi_hat) <= 2 * H^2 * (1 - exp(-lambda))")
    print("\nFor the given values:")
    print(f"H = {H}")
    print(f"lambda = {lmbda}")
    
    # Print the equation with numbers substituted
    print("\nThe calculation is:")
    term1 = 2
    term2 = H**2
    term3_val = 1 - math.exp(-lmbda)
    
    print(f"{term1} * {H}^2 * (1 - exp(-{lmbda})) = {term1} * {term2} * {term3_val:.4f} = {bound:.4f}")


# --- User-defined variables ---
# Please change these values to match your specific problem.
H_val = 20      # Example value for horizon H
lambda_val = 0.1  # Example value for hyperparameter lambda

# Calculate and print the bound
calculate_performance_bound(H_val, lambda_val)