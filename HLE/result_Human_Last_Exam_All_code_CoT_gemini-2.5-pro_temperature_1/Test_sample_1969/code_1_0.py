def display_shapley_formula():
    """
    Displays the derived formula for the Shapley value c_k for player p_k.
    The formula is presented term by term to show its structure.
    """

    print("The exact amount of money c_k for player p_k is given by the formula:")
    
    # Define the components of the formula
    sum_i = "n*(n+1)/2"
    sum_i_sq = "n*(n+1)*(2*n+1)/6"

    # Construct the terms of the formula
    term1 = f"k * ({sum_i})^3"
    term2 = f"k^2 * ({sum_i})^2"
    term3 = f"k * ({sum_i}) * ({sum_i_sq})"
    
    # Print the full formula
    print(f"\nc_k = {term1} - {term2} + {term3}\n")
    
    print("Here is the final equation with each number shown explicitly:")
    
    # Print the first term with its numbers
    print("Term 1: k * (n * (n + 1) / 2)^3")
    
    # Print the second term with its numbers
    print("Term 2: -1 * k^2 * (n * (n + 1) / 2)^2")
    
    # Print the third term with its numbers
    print("Term 3: k * (n * (n + 1) / 2) * (n * (n + 1) * (2*n + 1) / 6)")

display_shapley_formula()