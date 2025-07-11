def print_zero_count_formula():
    """
    Prints the formula for the least number of zeros of a vector field
    on a compact manifold M with a non-empty boundary ∂M.
    """
    
    # Define string representations for the mathematical terms
    chi_M = "χ(M)"
    chi_partial_M = "χ(∂M)"
    
    # The number in the denominator of the formula
    divisor = 2
    
    print("The least number of zeros a vector field on a compact manifold M with non-empty boundary ∂M can have is given by the formula:")
    print("") # Newline for readability
    
    # Print the formula in a human-readable format.
    # The f-string constructs the output string by embedding the variables.
    print(f"| {chi_M} - ( {chi_partial_M} / {divisor} ) |")
    print("")
    
    print(f"Here, '{chi_M}' is the Euler characteristic of the manifold M, and '{chi_partial_M}' is the Euler characteristic of its boundary ∂M.")
    print(f"The number '{divisor}' is printed as requested by the prompt.")

# Execute the function to print the formula
print_zero_count_formula()