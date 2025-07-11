def solve_magnetization_curve():
    """
    This function prints the analytical expression for the initial magnetization
    curve of a superconducting bar in the critical-state model.
    """
    
    # Define the variables for clarity in the printed expression
    M = "M"   # Magnetization
    H = "H"   # Applied magnetic field
    Jc = "Jc" # Constant critical current density
    b = "b"   # Half-thickness of the bar
    a = "a"   # Half-width of the bar
    
    # Explain the result and assumptions
    print("Based on the critical-state model for a thin slab in a parallel field,")
    print("the analytical expression for the initial magnetization M as a function of the")
    print("applied field H is derived.")
    print("-" * 60)
    
    # Construct and print the final analytical expression
    # We output each term of the equation as requested.
    print("Final Equation:")
    print(f"  {M}({H}) = - {H} * (1 - {H} / (2 * {Jc} * {b}))")
    print("-" * 60)
    
    # Explain the range of validity and the meaning of the variables
    print("This expression is valid for an increasing field H in the range:")
    print(f"  0 <= {H} <= {Jc} * {b}")
    print("\nWhere:")
    print(f"  {M}: The magnetization of the superconductor.")
    print(f"  {H}: The strength of the uniform applied magnetic field.")
    print(f"  {Jc}: The constant critical current density of the material.")
    print(f"  {b}: The half-thickness of the rectangular cross-section.")
    print(f"\nThe half-width '{a}' does not appear in the final expression because the")
    print(f"simplifying condition 'b << a' allows us to model the field penetration")
    print(f"as a one-dimensional problem across the smallest dimension (the thickness '2b').")

solve_magnetization_curve()

# The final answer in the required format
final_expression = "M(H) = -H * (1 - H / (2 * Jc * b))"
# Because the user wants an expression, not a single value, we'll format the output this way.
print(f"\n<<<M(H) = -H * (1 - H / (2 * Jc * b))>>>")