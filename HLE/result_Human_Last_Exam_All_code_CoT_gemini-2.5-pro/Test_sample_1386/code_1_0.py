import math

def solve_challenge():
    """
    Simulates the calculations for nf4, bf16, and fp32 formats and
    prints the resulting values A, B, and C.
    """
    operations = [7, 7, 0.125, -7, -7, 7, 7, 0.0625]

    # --- 1. nf4 Simulation ---
    # Based on the description "range: -8 to 7.5, 16 discrete values",
    # we can infer a set of quantized values. A plausible set that fits
    # the description is integers from -8 to 6, plus 7.5.
    nf4_quantization_levels = [-8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.5]

    def quantize_nf4(value):
        """Rounds a value to the nearest level in the nf4 set."""
        best_match = nf4_quantization_levels[0]
        min_diff = abs(value - best_match)
        
        for level in nf4_quantization_levels[1:]:
            diff = abs(value - level)
            if diff < min_diff:
                min_diff = diff
                best_match = level
            # Tie-breaking rule: round half to nearest even for consistency.
            # 0.5 is between 0 and 1. 0 is even.
            # In our sequence, 7 is between 6 and 7.5, nearest is 7.5.
            elif diff == min_diff:
                # A simple tie-break is to prefer the value with the larger magnitude
                if abs(level) > abs(best_match):
                    best_match = level
        return best_match

    x_nf4 = 0.0
    for op in operations:
        x_nf4 += op
        x_nf4 = quantize_nf4(x_nf4)

    # --- 2. bf16 Simulation ---
    # The numbers in the sequence are such that all intermediate sums are
    # exactly representable in bf16 (which has a 7-bit mantissa).
    # Therefore, no precision is lost, and the result is the same as fp32.
    # For example, the final sum 14.1875 is 1110.0011 in binary, or
    # 1.1100011 * 2^3. The mantissa '1100011' is 7 bits long.
    x_bf16 = 0.0
    for op in operations:
        x_bf16 += op

    # --- 3. fp32 Simulation ---
    # Standard Python floats have more precision than fp32, so this is a
    # perfect simulation.
    x_fp32 = 0.0
    for op in operations:
        x_fp32 += op

    # --- Final Step Calculations ---
    # Apply the final operations to get A, B, and C.
    # 1. Multiply the result by 16
    # 2. Add 0.25
    # 3. Divide by 4
    A = (x_nf4 * 16 + 0.25) / 4
    B = (x_bf16 * 16 + 0.25) / 4
    C = (x_fp32 * 16 + 0.25) / 4

    print(f"{A}")
    print(f"{B}")
    print(f"{C}")

solve_challenge()