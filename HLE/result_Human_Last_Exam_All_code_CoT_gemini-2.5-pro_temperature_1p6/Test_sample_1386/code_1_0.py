import struct
import math

def simulate_gpu_formats():
    """
    Simulates a sequence of arithmetic operations on three different
    GPU number formats: nf4, bf16, and fp32.
    """

    # The sequence of numbers to add
    sequence = [7, 7, 0.125, -7, -7, 7, 7, 0.0625]

    # --- 1. nf4 Simulation (Value A) ---
    # nf4 is simulated by clipping the result to its specified range [-8, 7.5]
    # after each operation. This is its most dominant characteristic.
    def quantize_nf4(value):
        return max(-8.0, min(7.5, value))

    x_nf4 = 0.0
    for num in sequence:
        x_nf4 = quantize_nf4(x_nf4 + num)

    # --- 2. bf16 Simulation (Value B) ---
    # bf16 is simulated by truncating a 32-bit float's mantissa to 7 bits
    # after each operation, which mimics the lower precision.
    def to_bfloat16_val(val):
        # Pack float to its 32-bit integer representation
        f_bits = struct.unpack('>I', struct.pack('>f', float(val)))[0]
        # Clear the lower 16 bits of the mantissa to simulate bfloat16
        f_bits &= 0xFFFF0000
        # Unpack back to a float
        return struct.unpack('>f', struct.pack('>I', f_bits))[0]

    x_bf16 = 0.0
    for num in sequence:
        x_bf16 = to_bfloat16_val(x_bf16 + num)

    # --- 3. fp32 Simulation (Value C) ---
    # fp32 is the standard Python float. No special handling needed.
    x_fp32 = 0.0
    for num in sequence:
        x_fp32 += num

    # --- Final operations for A, B, and C ---
    # These operations are applied to the result from each simulation.
    def final_transform(value):
        return (value * 16 + 0.25) / 4

    A = final_transform(x_nf4)
    B = final_transform(x_bf16)
    C = final_transform(x_fp32)
    
    print("The final values for A, B, and C are:")
    print(f"A (from nf4) = {A}")
    print(f"B (from bf16) = {B}")
    print(f"C (from fp32) = {C}")
    print("\nFinal equation to evaluate: ceil((B - C - A) * 10)")
    print(f"Plugging in the values: ceil(({B} - {C} - {A}) * 10)")

simulate_gpu_formats()