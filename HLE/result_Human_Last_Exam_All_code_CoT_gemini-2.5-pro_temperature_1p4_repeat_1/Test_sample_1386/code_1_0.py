import math

def solve_gpu_challenge():
    """
    Solves the GPU Number Format Detection Challenge by simulating the calculation
    in nf4, bf16, and fp32 formats.
    """
    sequence = [7, 7, 0.125, -7, -7, 7, 7, 0.0625]

    # --- 1. nf4 Simulation (A) ---
    # The nf4 format is defined by its small range [-8, 7.5].
    # The primary source of error will be clamping to this range.
    x_nf4 = 0.0
    nf4_min, nf4_max = -8.0, 7.5
    for num in sequence:
        x_nf4 += num
        # Clamp the result to the nf4 range
        if x_nf4 > nf4_max:
            x_nf4 = nf4_max
        elif x_nf4 < nf4_min:
            x_nf4 = nf4_min
    
    # Final operations for A
    res_nf4 = x_nf4
    A = (res_nf4 * 16 + 0.25) / 4

    # --- 2. bf16 Simulation (B) ---
    # The bf16 format has a 7-bit mantissa, leading to rounding errors.
    def round_to_bf16(val):
        """Simulates rounding to bf16 precision."""
        if val == 0:
            return 0.0
        # ULP = 2^(floor(log2(|val|)) - 7)
        exponent = math.floor(math.log2(abs(val)))
        mantissa_bits = 7
        ulp = 2**(exponent - mantissa_bits)
        return round(val / ulp) * ulp

    x_bf16 = 0.0
    for num in sequence:
        x_bf16 = round_to_bf16(x_bf16 + num)
    
    # Final operations for B, with rounding at each step
    res_bf16 = x_bf16
    val_b_1 = round_to_bf16(res_bf16 * 16)
    val_b_2 = round_to_bf16(val_b_1 + 0.25)
    B = round_to_bf16(val_b_2 / 4)

    # --- 3. fp32 Simulation (C) ---
    # Python floats (fp64) have more than enough precision to simulate fp32 accurately.
    x_fp32 = 0.0
    for num in sequence:
        x_fp32 += num
    
    # Final operations for C
    res_fp32 = x_fp32
    C = (res_fp32 * 16 + 0.25) / 4

    # --- Print Results ---
    print("Simulating the calculations for each format...\n")
    print(f"Value A (from nf4 simulation): {A}")
    print(f"Value B (from bf16 simulation): {B}")
    print(f"Value C (from fp32 simulation): {C}")
    print("\nAssembling the final expression for mental calculation:")
    print(f"ceil((B - C - A) * 10)")
    print(f"ceil(({B} - {C} - {A}) * 10)")

solve_gpu_challenge()