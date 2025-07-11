def print_final_answer():
    """
    This function prints the final derived expressions for the magnetic field H,
    which corresponds to choice E.
    """
    
    # Header
    print("The final expressions for the magnetic field H(r, θ) are given by:")
    
    # Define the equations as multiline strings for clarity
    h_field_case = "H(r, θ) = "
    
    h_field_inside = r"""
    { 
    {  (2 * μ₀ / μ) * (K₀ / (1 + (2 * μ₀ / μ))) * ẑ      for 0 < r < R
    {"""

    h_field_outside = r"""
    {  (K₀ / (1 + (2 * μ₀ / μ))) * (R³ / r³) * (2*cos(θ)*r̂ + sin(θ)*θ̂)      for R < r < ∞
    {"""

    # Printing the formatted answer
    print(h_field_case)
    # The string replacements are to align the braces properly for a combined view
    print(h_field_inside.replace("{", " ", 1))
    print(h_field_outside.replace("{", " ", 1).rstrip())


# Execute the function to print the solution
print_final_answer()