def explain_theoretical_stress():
    """
    This function explains the calculation of the theoretical stress at the tip of a sharp wedge.
    """
    print("To find the theoretical stress at the tip of the wedge (Point A), we use the principles of stress concentration.")
    print("\nThe governing equation for the maximum stress at a notch is:")
    # We represent the equation using strings.
    sigma_max = "σ_max"
    sigma_y = "σ_y"
    a = "a"
    rho = "ρ"
    
    # Print the equation with its components
    print(f"\n  {sigma_max} = {sigma_y} * (1 + 2 * ({a} / {rho}))\n")
    
    print("Where:")
    print(f"  {sigma_max} is the maximum stress at the tip of the wedge.")
    print(f"  {sigma_y} is the applied uniaxial stress.")
    print(f"  {a} is the depth of the notch.")
    print(f"  {rho} is the radius of curvature of the notch tip.")
    
    print("\nFor a perfectly sharp theoretical wedge, the tip is a point, so its radius of curvature 'ρ' is 0.")
    
    print("\nAnalyzing the equation as 'ρ' approaches 0:")
    print(f"  The term '({a} / {rho})' approaches infinity.")
    print(f"  Therefore, the entire expression for '{sigma_max}' also approaches infinity.")
    
    print("\nConclusion: The theoretical stress at the tip of a perfectly sharp wedge is infinite.")
    print(f"\nFinal Equation Result: {sigma_max} = ∞")

explain_theoretical_stress()