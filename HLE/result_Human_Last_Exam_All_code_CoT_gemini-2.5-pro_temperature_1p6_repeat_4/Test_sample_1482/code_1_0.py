import math

def solve_it3_vertical_cross_section():
    """
    This function derives and presents the mathematical formulation for the
    vertical cross-section of an Interval Type-3 Membership Function (IT3 MF)
    based on a Gaussian paradigm.
    """

    # --- 1. General Formula ---
    # Define symbols using unicode for clarity
    c = "c"
    sigma_1 = "σ₁"
    sigma_2 = "σ₂"
    x_sym = "x"
    u_sym = "u"
    mu_A_x_u = f"μ_A({x_sym})({u_sym})"

    # The formula describes the membership grade (alpha) of a primary membership
    # value 'u' for a given input 'x'.
    print("### Mathematical Formulation ###\n")
    print("The vertical cross-section of an upper Interval Type-3 Gaussian membership function,")
    print("where the uncertainty is in the standard deviation σ ∈ [σ₁, σ₂], is given by:")
    
    # Print the general formula
    formula_str = f"{mu_A_x_u} = (1 / ({sigma_2} - {sigma_1})) * [ |{x_sym} - {c}| / sqrt(-2 * ln({u_sym})) - {sigma_1} ]"
    print("\n" + formula_str + "\n")
    print("-" * 30)

    # --- 2. Example Calculation ---
    print("\n### Example Calculation ###\n")
    # Define concrete values for the parameters
    c_val = 5.0
    sigma_1_val = 1.0
    sigma_2_val = 2.0
    x_val = 6.0
    
    # The value for u must be within the range defined by σ₁ and σ₂ for the given x
    u_max = math.exp(-0.5 * ((x_val - c_val) / sigma_1_val)**2)
    u_min = math.exp(-0.5 * ((x_val - c_val) / sigma_2_val)**2)
    
    # We choose a u value within this valid range [u_min, u_max]
    # For this example, let's pick a value in the middle.
    u_val = (u_max + u_min) / 2
    
    print(f"Given the following parameter values:")
    print(f"  Center {c} = {c_val}")
    print(f"  Standard Deviation Interval [{sigma_1}, {sigma_2}] = [{sigma_1_val}, {sigma_2_val}]")
    print(f"  Primary Input Variable {x_sym} = {x_val}")
    print(f"  Primary Membership Value {u_sym} = {u_val:.4f}\n")
    
    # --- 3. Substitute numbers into the formula ---
    print("Substituting these values into the equation:\n")
    
    # Print the equation with numbers
    # Each number in the final equation is outputted explicitly.
    numerator_abs_val = f"|{x_val} - {c_val}|"
    denominator_sqrt = f"sqrt(-2 * ln({u_val:.4f}))"
    denominator_diff = f"({sigma_2_val} - {sigma_1_val})"

    equation_with_numbers = f"μ_A({x_val})({u_val:.4f}) = (1 / {denominator_diff}) * [ {numerator_abs_val} / {denominator_sqrt} - {sigma_1_val} ]"
    print(equation_with_numbers)

    # --- 4. Calculate the final result ---
    try:
        # Step-by-step calculation
        abs_diff = abs(x_val - c_val)
        ln_u = math.log(u_val)
        sqrt_term = math.sqrt(-2 * ln_u)
        
        # This intermediate value is sigma
        sigma_val = abs_diff / sqrt_term
        
        # This is the final result (alpha)
        result = (1 / (sigma_2_val - sigma_1_val)) * (sigma_val - sigma_1_val)

        print("\nIntermediate calculation step (solving for σ):")
        print(f"σ = {abs_diff:.1f} / sqrt(-2 * {ln_u:.4f}) = {sigma_val:.4f}")

        print("\nFinal calculation:")
        print(f"μ_A({x_val})({u_val:.4f}) = (1 / {sigma_2_val - sigma_1_val:.1f}) * [ {sigma_val:.4f} - {sigma_1_val:.1f} ] = {result:.4f}")

    except ValueError as e:
        print(f"\nCalculation Error: {e}. This typically occurs if 'u' is outside the valid range for the given 'x'.")

# Run the function
solve_it3_vertical_cross_section()