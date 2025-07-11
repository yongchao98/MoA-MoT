import sys

def calculate_alpha_scaling():
    """
    Calculates and prints the scaling relation for the specific heat
    critical exponent alpha with spatial dimensionality d below d=4,
    for a scalar field theory (n=1).
    """

    # In the context of the O(n) model, a scalar field theory corresponds to n=1.
    n = 1

    # The upper critical dimension for this model is d_c = 4.
    d_c = 4

    # The epsilon-expansion provides a formula for alpha to first order in epsilon = d_c - d.
    # The general formula is: α ≈ ((4 - n) / (n + 8)) * ε

    # We calculate the numerator and denominator of the coefficient for n=1.
    numerator = 4 - n
    denominator = n + 8

    # Calculate the coefficient as a floating-point number.
    coefficient = numerator / denominator

    # --- Outputting the Derivation and Final Equation ---

    print("The quantitative scaling of the specific heat exponent (α) with dimensionality (d) is found via the ε-expansion, where ε = 4 - d.")
    print("\nFor a scalar field (O(n) model with n=1), the first-order formula is:")
    print("α ≈ ((4 - n) / (n + 8)) * ε")

    print("\nSubstituting the numbers for this model:")
    print(f"n = {n}")
    print(f"d_c = {d_c} (upper critical dimension)")

    print("\nThe scaling relation with all numbers explicitly shown is:")
    # The user requested that each number in the final equation be output.
    print(f"α ≈ (({d_c} - {n}) / ({n} + {denominator})) * ({d_c} - d)")
    print(f"α ≈ ({numerator} / {denominator}) * ({d_c} - d)")
    
    print("\nThis simplifies to the final quantitative relationship:")
    print(f"α ≈ {coefficient:.4f} * ({d_c} - d)")


# Execute the function
calculate_alpha_scaling()

# --- Final Answer Formatting ---
# The final answer is the derived formula.
# We will construct it as a string to be embedded in the final output format.
final_answer_formula = f"α ≈ (1/3) * (4 - d)"

# Appending the final answer in the specified format to the output.
# We add newlines for clean separation.
sys.stdout.write('\n' * 2)
sys.stdout.write(f'<<<{final_answer_formula}>>>')