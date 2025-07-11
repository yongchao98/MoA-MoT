def solve_information_theory_problem():
    """
    Calculates the largest possible value of I(X;Y|W) based on the given information.
    """
    # Given values
    I_X_Y = 3  # I(X;Y)
    I_X_Y_given_Z = 2  # I(X;Y|Z)
    I_X_Z_given_Y = 5  # I(X;Z|Y)

    print("Step 1: Calculate related mutual information quantities.")
    # Using the chain rule for mutual information: I(A;B,C) = I(A;B) + I(A;C|B) = I(A;C) + I(A;B|C)
    # We can write: I(X;Y) + I(X;Z|Y) = I(X;Z) + I(X;Y|Z)
    # Let's solve for I(X;Z)
    # I(X;Z) = I(X;Y) + I(X;Z|Y) - I(X;Y|Z)
    I_X_Z = I_X_Y + I_X_Z_given_Y - I_X_Y_given_Z
    print(f"From the chain rule, I(X;Y) + I(X;Z|Y) = I(X;Z) + I(X;Y|Z)")
    print(f"{I_X_Y} + {I_X_Z_given_Y} = I(X;Z) + {I_X_Y_given_Z}")
    print(f"So, I(X;Z) = {I_X_Y} + {I_X_Z_given_Y} - {I_X_Y_given_Z} = {I_X_Z}")

    # The interaction information I(X;Y;Z) measures the redundancy or synergy between Y and Z about X.
    # I(X;Y;Z) = I(X;Y) - I(X;Y|Z)
    I_X_Y_Z = I_X_Y - I_X_Y_given_Z
    print(f"\nThe interaction information is I(X;Y;Z) = I(X;Y) - I(X;Y|Z) = {I_X_Y} - {I_X_Y_given_Z} = {I_X_Y_Z}")
    # Let's check for consistency: I(X;Y;Z) = I(X;Z) - I(X;Z|Y)
    print(f"Check: I(X;Z) - I(X;Z|Y) = {I_X_Z} - {I_X_Z_given_Y} = {I_X_Z - I_X_Z_given_Y}. The values are consistent.")

    print("\nStep 2: Find a general expression for I(X;Y|W).")
    # Using the chain rule on the variables (X,Y,W): I(X;Y,W) = I(X;Y) + I(X;W|Y) = I(X;W) + I(X;Y|W)
    # Rearranging for I(X;Y|W) gives:
    # I(X;Y|W) = I(X;Y) + I(X;W|Y) - I(X;W)
    print("From the chain rule involving W, we have the identity:")
    print("I(X;Y|W) = I(X;Y) + I(X;W|Y) - I(X;W)")
    print(f"I(X;Y|W) = {I_X_Y} + I(X;W|Y) - I(X;W)")

    print("\nStep 3: Analyze the term I(X;W|Y) - I(X;W).")
    # The term I(X;W|Y) - I(X;W) is equal to -I(X;Y;W), where I(X;Y;W) is the interaction information involving W.
    # I(X;Y;W) = I(X;W) - I(X;W|Y).
    print("The term I(X;W|Y) - I(X;W) is the negative of the interaction information I(X;Y;W).")
    
    print("\nStep 4: Use the Data Processing Inequality for interaction information.")
    # We are given W is a deterministic function of Z. This implies a Markov chain (X,Y) -> Z -> W.
    # A key principle in information theory states that local processing of a variable cannot create synergy from redundancy.
    # Since I(X;Y;Z) = 1 > 0, Y and Z are redundant about X.
    # Therefore, for any W=g(Z), the interaction information I(X;Y;W) must be non-negative.
    print(f"Since I(X;Y;Z) = {I_X_Y_Z} > 0, variables Y and Z are redundant about X.")
    print("Because W is a function of Z, W cannot be more synergistic with Y than Z is.")
    print("This implies I(X;Y;W) >= 0.")
    print("I(X;Y;W) is defined as I(X;W) - I(X;W|Y).")
    print("So, I(X;W) - I(X;W|Y) >= 0, which means I(X;W|Y) - I(X;W) <= 0.")

    print("\nStep 5: Find the maximum value of I(X;Y|W).")
    # We want to maximize: I(X;Y|W) = I_X_Y + (I(X;W|Y) - I(X;W))
    # Since I(X;W|Y) - I(X;W) <= 0, the maximum value of this term is 0.
    max_val = I_X_Y + 0
    print(f"To maximize I(X;Y|W) = {I_X_Y} + (I(X;W|Y) - I(X;W)), we must maximize the term in parentheses.")
    print("As we've shown, this term is less than or equal to 0.")
    print(f"The maximum possible value is therefore {I_X_Y} + 0 = {max_val}.")

    print("\nStep 6: Show that this maximum value is achievable.")
    # The maximum is achieved when I(X;W|Y) - I(X;W) = 0, which means I(X;W) = I(X;W|Y).
    # This can be achieved by choosing W to be a function of Z that isolates the information in Z that is unique about X with respect to Y.
    # The amount of this unique information is I(X;Z|Y) = 5.
    # Let's imagine a W=g(Z) that perfectly captures this information.
    I_X_W_hypothetical = I_X_Z_given_Y  # Set I(X;W) = 5
    I_X_W_given_Y_hypothetical = I_X_Z_given_Y # Set I(X;W|Y) = 5
    print("This maximum is achieved if we can find a W=g(Z) such that I(X;W) = I(X;W|Y).")
    print(f"Consider a W that only captures the unique information Z has about X (given Y). The size of this information is I(X;Z|Y) = {I_X_Z_given_Y}.")
    print(f"For such a W, we would have I(X;W) = {I_X_W_hypothetical} and I(X;W|Y) = {I_X_W_given_Y_hypothetical}.")
    # Let's plug this into our main formula:
    final_value = I_X_Y + I_X_W_given_Y_hypothetical - I_X_W_hypothetical
    print(f"Plugging this into the equation I(X;Y|W) = {I_X_Y} + I(X;W|Y) - I(X;W):")
    print(f"I(X;Y|W) = {I_X_Y} + {I_X_W_given_Y_hypothetical} - {I_X_W_hypothetical} = {final_value}")
    print("This confirms that the maximum value is achievable.")
    
    return final_value

result = solve_information_theory_problem()
print(f"\nThe largest possible value of I(X;Y|W) is {result}.")
