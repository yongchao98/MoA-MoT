import struct
import math

# --- Helper functions to simulate number formats ---

def to_nf4(value):
    """
    Simulates the nf4 format by clamping the value to its specified range.
    This model assumes clamping is the dominant effect for the given sequence.
    """
    NF4_MIN = -8.0
    NF4_MAX = 7.5
    if value > NF4_MAX:
        return NF4_MAX
    elif value < NF4_MIN:
        return NF4_MIN
    return value

def to_bfloat16(value):
    """
    Simulates bfloat16 (bf16) rounding by converting a Python float (fp64)
    to fp32 and then truncating the mantissa. This captures the precision loss effect.
    """
    # Pack the float to a 32-bit integer representation
    try:
        bits = struct.unpack('<I', struct.pack('<f', float(value)))[0]
    except OverflowError:
        return float('inf') if value > 0 else float('-inf')
    
    # Round to nearest by adding 0x8000, then mask the lower 16 bits to simulate truncation
    # This effectively truncates the 23-bit fp32 mantissa to bf16's 7-bit mantissa region
    bits = bits + 0x8000
    bits = bits & 0xFFFF0000

    # Unpack the bits back to a float
    return struct.unpack('<f', struct.pack('<I', bits))[0]

# --- Main calculation ---

# The sequence of numbers to be added
add_sequence = [7, 7, 0.125, -7, -7, 7, 7, 0.0625]

# Case 1: nf4 calculation (Value A)
val_a = 0.0
for num in add_sequence:
    val_a = to_nf4(val_a + num)

# Apply final operations for nf4
val_a = to_nf4(val_a * 16)
val_a = to_nf4(val_a + 0.25)
val_a = to_nf4(val_a / 4)

# Case 2: bf16 calculation (Value B)
val_b = 0.0
for num in add_sequence:
    val_b = to_bfloat16(val_b + num)

# Apply final operations for bf16
val_b = to_bfloat16(val_b * 16)
val_b = to_bfloat16(val_b + 0.25)
val_b = to_bfloat16(val_b / 4)

# Case 3: fp32 calculation (Value C)
val_c = 0.0
for num in add_sequence:
    val_c += num

# Apply final operations for fp32
val_c = val_c * 16
val_c = val_c + 0.25
val_c = val_c / 4

# --- Output the results ---
print("The final values for the equation ceil((B - C - A) * 10) are:")
print(f"A (nf4) = {val_a}")
print(f"B (bf16) = {val_b}")
print(f"C (fp32) = {val_c}")

# The user is asked to perform the final calculation mentally.
# B - C - A = 56.75 - 56.8125 - 1.875 = -1.9375
# (B - C - A) * 10 = -19.375
# ceil(-19.375) = -19