def calculate_minimax_risk(n):
    """
    Calculates the minimax risk for estimating the parameter theta of a
    Binomial distribution based on n i.i.d. observations of a Bin(n, theta)
    random variable.

    Args:
        n (int): The number of observations (and the number of trials in each
                 observation). Must be a positive integer.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return

    # Based on the problem description, the sufficient statistic follows a
    # Bin(n^2, theta) distribution. Let N = n^2.
    # The general formula for the minimax risk is 1 / (4 * (sqrt(N) + 1)^2).
    # For N = n^2, this simplifies to 1 / (4 * (n + 1)^2).

    print(f"For n = {n}:")
    print(f"The problem is equivalent to estimating theta from a single observation T ~ Bin(N, theta), where N = n*n = {n**2}.")
    print("The minimax risk formula simplifies to: 1 / (4 * (n + 1)^2)")
    print("-" * 30)

    # Step-by-step calculation based on the simplified formula
    print("Calculation breakdown:")
    
    # 1. Calculate n + 1
    val_n_plus_1 = n + 1
    print(f"Step 1: Calculate (n + 1) -> {n} + 1 = {val_n_plus_1}")

    # 2. Calculate (n + 1)^2
    val_n_plus_1_squared = val_n_plus_1**2
    print(f"Step 2: Calculate (n + 1)^2 -> ({val_n_plus_1})^2 = {val_n_plus_1_squared}")

    # 3. Calculate the denominator 4 * (n + 1)^2
    denominator = 4 * val_n_plus_1_squared
    print(f"Step 3: Calculate 4 * (n + 1)^2 -> 4 * {val_n_plus_1_squared} = {denominator}")
    print("-" * 30)
    
    # 4. Calculate the final risk
    minimax_risk = 1 / denominator
    print("Final minimax risk:")
    print(f"Risk = 1 / {denominator} = {minimax_risk}")


# --- User Execution ---
# You can change the value of n to see the result for a different case.
n_value = 5
calculate_minimax_risk(n_value)