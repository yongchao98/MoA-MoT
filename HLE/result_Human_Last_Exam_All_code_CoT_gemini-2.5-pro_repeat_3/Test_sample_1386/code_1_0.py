import math

def simulate_gpu_formats():
    """
    Simulates a sequence of floating-point operations on three different
    GPU number formats (nf4, bf16, fp32) and calculates the final values A, B, and C.
    """
    numbers = [7, 7, 0.125, -7, -7, 7, 7, 0.0625]

    # --- nf4 Simulation ---
    # The key feature of nf4 here is the limited range [-8, 7.5].
    # Any value outside this range is clamped.
    val_nf4 = 0.0
    nf4_min = -8.0
    nf4_max = 7.5
    for num in numbers:
        val_nf4 += num
        if val_nf4 > nf4_max:
            val_nf4 = nf4_max
        elif val_nf4 < nf4_min:
            val_nf4 = nf4_min
    res_nf4 = val_nf4

    # --- bf16 and fp32 Simulation ---
    # For this specific sequence, neither bf16 (7-bit mantissa) nor fp32 (23-bit mantissa)
    # will lose precision. The intermediate and final sums are exactly representable in both.
    # Therefore, their results will be identical to a standard high-precision calculation.
    val_high_precision = 0.0
    for num in numbers:
        val_high_precision += num
    res_bf16 = val_high_precision
    res_fp32 = val_high_precision

    # --- Final Calculations ---
    # 1. Multiply the result by 16
    # 2. Add 0.25
    # 3. Divide by 4
    A = (res_nf4 * 16 + 0.25) / 4
    B = (res_bf16 * 16 + 0.25) / 4
    C = (res_fp32 * 16 + 0.25) / 4

    print("Intermediate result from addition sequence:")
    print(f"nf4: {res_nf4}")
    print(f"bf16: {res_bf16}")
    print(f"fp32: {res_fp32}\n")

    print("Final Values:")
    print(f"A = (({res_nf4} * 16) + 0.25) / 4 = {A}")
    print(f"B = (({res_bf16} * 16) + 0.25) / 4 = {B}")
    print(f"C = (({res_fp32} * 16) + 0.25) / 4 = {C}")

    # The final step is performed mentally by the user based on the output.
    # The value is ceil((B - C - A) * 10)
    # B - C - A = 56.8125 - 56.8125 - 30.0625 = -30.0625
    # (B - C - A) * 10 = -300.625
    # ceil(-300.625) = -300

simulate_gpu_formats()