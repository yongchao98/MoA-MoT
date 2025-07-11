import sympy

def calculate_cost_ratio():
    """
    Calculates and prints the symbolic formula for the ratio of additional
    computational cost.
    """
    # Define the symbolic variables used in the formula
    L, d_model, d_attn, d_ff = sympy.symbols('L d_model d_attn d_ff')

    # The numerator of the derived ratio is 3 * L * d_attn
    numerator = 3 * L * d_attn

    # The denominator of the derived ratio is
    # 2 * d_model * (2 * d_attn + d_ff) + L * d_attn
    denominator = 2 * d_model * (2 * d_attn + d_ff) + L * d_attn

    # The problem asks to output each number in the final equation.
    # Since the answer is a symbolic formula, we will print the formula's structure.
    print("The best estimate for the ratio of additional computational cost is:")
    print("      3 * L * d_attn")
    print("---------------------------------------------")
    print("2 * d_model * (2 * d_attn + d_ff) + L * d_attn")

calculate_cost_ratio()