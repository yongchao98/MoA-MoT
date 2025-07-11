import struct
import math

# --- Helper Functions for Format Simulation ---

def simulate_bf16_conversion(float_val: float) -> float:
    """
    Simulates float32 to bfloat16 conversion using bitwise truncation.
    This round-towards-zero method is sufficient for the numbers in this problem.
    """
    # Pack the float into 4 bytes (32 bits)
    packed_float = struct.pack('<f', float_val)
    # Unpack as an unsigned integer
    int_val = struct.unpack('<I', packed_float)[0]
    
    # Truncate the last 16 bits (the low-precision part of fp32 mantissa)
    int_val = int_val >> 16
    
    # Shift back to restore exponent and sign position
    int_val = int_val << 16
    
    # Repack as bytes and unpack as a float
    packed_val = struct.pack('<I', int_val)
    return struct.unpack('<f', packed_val)[0]

def simulate_nf4_quantization(float_val: float) -> float:
    """
    Simulates the nf4 format, focusing on its dominant feature for this problem:
    range clipping.
    """
    NF4_MIN = -8.0
    NF4_MAX = 7.5
    
    # Clamp the value to the specified nf4 range
    if float_val > NF4_MAX:
        return NF4_MAX
    if float_val < NF4_MIN:
        return NF4_MIN
    return float_val

# --- Main Calculation Logic ---

# The sequence of numbers to be added
sequence = [7, 7, 0.125, -7, -7, 7, 7, 0.0625]

# 1. Calculation for fp32 (Value C)
# Standard Python floats are fp32. No special handling needed.
c = 0.0
for num in sequence:
    c += num
c = (c * 16 + 0.25) / 4

# 2. Calculation for bf16 (Value B)
# Each arithmetic operation's result is converted to bf16.
b = 0.0
for num in sequence:
    b = simulate_bf16_conversion(b + num)
    
# The final operations are also subject to bf16 precision.
b = simulate_bf16_conversion(b * 16)
b = simulate_bf16_conversion(b + 0.25)
b = simulate_bf16_conversion(b / 4)

# 3. Calculation for nf4 (Value A)
# Each addition is subject to nf4's range limitations.
a = 0.0
for num in sequence:
    a = simulate_nf4_quantization(a + num)
    
# The final operations are performed in a higher precision.
a = (a * 16 + 0.25) / 4


# --- Output the Results ---
print("GPU Number Format Detection Results:")
print(f"A (nf4 final value): {a}")
print(f"B (bf16 final value): {b}")
print(f"C (fp32 final value): {c}")
print("\nFinal Equation:")
print("ceil((B - C - A) * 10)")
print("Breaking down the expression:")
expression_result = (b - c - a) * 10
print(f"B - C - A  = {b} - {c} - {a} = {b - c - a}")
print(f"(B - C - A) * 10 = {expression_result}")
print(f"ceil({expression_result}) = {math.ceil(expression_result)}")
