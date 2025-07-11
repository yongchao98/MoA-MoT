import struct
import math

def quantize_bfloat16(f_val):
    """
    Simulates a float32 to bfloat16 to float32 conversion.
    This is done by truncating the 16 least significant bits of the mantissa.
    """
    # Use struct to get the 32-bit integer representation of the float
    packed = struct.pack('>f', f_val)
    int_val = struct.unpack('>I', packed)[0]
    
    # Truncate the last 16 bits (the end of the mantissa)
    int_val &= 0xFFFF0000
    
    # Pack the integer back into a float
    packed_again = struct.pack('>I', int_val)
    return struct.unpack('>f', packed_again)[0]

def quantize_nf4(f_val):
    """
    Simulates nf4 quantization by clipping the value to the specified range.
    """
    NF4_MIN = -8.0
    NF4_MAX = 7.5
    return max(NF4_MIN, min(f_val, NF4_MAX))

def solve():
    """
    Runs the simulation for fp32, bf16, and nf4 formats.
    """
    sequence = [7, 7, 0.125, -7, -7, 7, 7, 0.0625]

    # --- C: fp32 Simulation ---
    c_val = 0.0
    for num in sequence:
        c_val += num
    c_val *= 16
    c_val += 0.25
    c_val /= 4
    
    # --- B: bf16 Simulation ---
    b_val = 0.0
    b_val = quantize_bfloat16(b_val)
    for num in sequence:
        b_val = quantize_bfloat16(b_val + num)
    b_val = quantize_bfloat16(b_val * 16)
    b_val = quantize_bfloat16(b_val + 0.25)
    b_val = quantize_bfloat16(b_val / 4)

    # --- A: nf4 Simulation ---
    a_val = 0.0
    a_val = quantize_nf4(a_val)
    for num in sequence:
        a_val = quantize_nf4(a_val + num)
    a_val = quantize_nf4(a_val * 16)
    a_val = quantize_nf4(a_val + 0.25)
    a_val = quantize_nf4(a_val / 4)

    print("Final values for the equation ceil((B-C-A)*10):")
    print(f"A (nf4) = {a_val}")
    print(f"B (bf16) = {b_val}")
    print(f"C (fp32) = {c_val}")
    
    print("\nThe equation with the computed values is:")
    print(f"ceil(({b_val} - {c_val} - {a_val}) * 10)")
    
solve()