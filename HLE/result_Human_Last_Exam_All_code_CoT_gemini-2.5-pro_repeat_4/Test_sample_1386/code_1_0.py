import numpy as np
import math

def solve_gpu_format_challenge():
    """
    This script simulates a sequence of arithmetic operations in three different
    GPU number formats (nf4, bf16, fp32) to determine the final result for each.
    """
    
    # The sequence of numbers to be added
    numbers = [7, 7, 0.125, -7, -7, 7, 7, 0.0625]

    # --- Case 1: nf4 (Value A) ---
    # nf4 is a 4-bit format with a specified range of [-8, 7.5].
    # This implies very low precision and clipping for out-of-range values.
    # We simulate this by clipping results to the range and rounding to the
    # nearest integer after each operation.
    
    print("--- Calculating A for nf4 ---")
    
    def quantize_nf4(v):
        # 1. Clip the value to the specified range of nf4
        clipped = max(-8.0, min(7.5, v))
        # 2. Round to the nearest integer to simulate low precision.
        # Python's round() rounds .5 to the nearest even integer.
        # We use a custom implementation for clarity.
        if clipped > 0:
            quantized = math.floor(clipped + 0.5)
        else:
            quantized = math.ceil(clipped - 0.5)
        return float(quantized)

    x_nf4 = 0.0
    for num in numbers:
        x_nf4 = quantize_nf4(x_nf4 + num)
        
    intermediate_A = x_nf4
    
    # The final operations are also performed in nf4
    res_A_1 = quantize_nf4(intermediate_A * 16)
    res_A_2 = quantize_nf4(res_A_1 + 0.25)
    A = quantize_nf4(res_A_2 / 4)
    
    print(f"Intermediate value x = {intermediate_A}")
    print(f"Final equation: Q(Q(Q({intermediate_A} * 16) + 0.25) / 4)")
    print(f"Result A (nf4) = {A}")
    print("-" * 25)

    # --- Case 2: bf16 (Value B) ---
    # bf16 has the same range as fp32 but only 7 bits of precision.
    # This can lead to loss of precision when adding numbers of different magnitudes.
    # We use numpy's bfloat16 for an accurate simulation.
    
    print("--- Calculating B for bf16 ---")
    
    x_bf16 = np.bfloat16(0)
    for num in numbers:
        x_bf16 += np.bfloat16(num)

    intermediate_B = x_bf16
    
    # The final operations are also performed in bf16
    res_B_1 = intermediate_B * np.bfloat16(16)
    res_B_2 = res_B_1 + np.bfloat16(0.25) # Precision is lost here
    B = res_B_2 / np.bfloat16(4)
    B = float(B) # Convert back to standard float for printing
    
    print(f"Intermediate value x = {intermediate_B}")
    print(f"Final equation: ({intermediate_B} * 16 + 0.25) / 4")
    print(f"Result B (bf16) = {B}")
    print("-" * 25)
    
    # --- Case 3: fp32 (Value C) ---
    # fp32 is the standard single-precision float. Python's floats are
    # double-precision (fp64), which is more than enough to get the correct fp32 result.
    
    print("--- Calculating C for fp32 ---")
    
    x_fp32 = 0.0
    for num in numbers:
        x_fp32 += num

    intermediate_C = x_fp32
    
    # The final operations
    C = (intermediate_C * 16 + 0.25) / 4
    
    print(f"Intermediate value x = {intermediate_C}")
    print(f"Final equation: ({intermediate_C} * 16 + 0.25) / 4")
    print(f"Result C (fp32) = {C}")
    print("-" * 25)

solve_gpu_format_challenge()