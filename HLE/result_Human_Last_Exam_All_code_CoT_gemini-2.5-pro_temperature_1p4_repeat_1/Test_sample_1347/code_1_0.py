def get_r0f_expression():
    """
    Constructs and prints the symbolic expression for R0f, the basic reproduction number for forest fires.
    """
    
    # Define the variables as strings for creating the formula
    b = "b"
    c = "c"
    pg = "pg"
    pt = "pt"
    gamma_t = "𝛾t"
    mu_t = "𝜇t"
    mu_g = "𝜇g"
    
    # Numerator of the R0f expression
    # Represents the combined rates and probabilities of transmission
    numerator = f"{b} * {c} * {pg} * {pt}"
    
    # Denominator of the R0f expression
    # Represents the product of the removal rates for burning trees and burning grass
    denominator = f"({gamma_t} + {mu_t}) * {mu_g}"
    
    # The complete expression for R0f
    expression = f"R0f = ({numerator}) / ({denominator})"
    
    # Print the final expression
    print(expression)

# Execute the function to see the output
get_r0f_expression()