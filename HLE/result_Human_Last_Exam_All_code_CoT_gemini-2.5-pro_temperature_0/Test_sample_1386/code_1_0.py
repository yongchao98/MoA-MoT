import struct
import math

def quantize_nf4(x):
    """
    Simulates a simplified nf4 format.
    Assumes 16 discrete values from -8 to 7 (4-bit signed integer).
    Values are rounded to the nearest integer and clamped to the [-8, 7] range.
    """
    return max(-8.0, min(7.0, round(x)))

def quantize_bf16(f):
    """
    Simulates bfloat16 quantization by truncating the mantissa of an fp32.
    It packs the float into bytes, truncates the last 16 bits, and unpacks it.
    """
    # Pack the float into its 32-bit big-endian representation
    b = struct.pack('>f', f)
    # Convert bytes to an integer, right-shift 16 bits to truncate, then shift back
    i = int.from_bytes(b, 'big')
    i = (i >> 16) << 16
    # Convert the modified integer back to bytes and then to a float
    b_new = i.to_bytes(4, 'big')
    return struct.unpack('>f', b_new)[0]

def solve():
    """
    Runs the calculation sequence for nf4, bf16, and fp32 formats.
    """
    sequence = [7, 7, 0.125, -7, -7, 7, 7, 0.0625]

    # --- nf4 Calculation (A) ---
    a = 0.0
    for num in sequence:
        a = quantize_nf4(a + num)
    a = quantize_nf4(a * 16)
    a = quantize_nf4(a + 0.25)
    a = quantize_nf4(a / 4)

    # --- bf16 Calculation (B) ---
    b = 0.0
    for num in sequence:
        b = quantize_bf16(b + num)
    b = quantize_bf16(b * 16)
    b = quantize_bf16(b + 0.25)
    b = quantize_bf16(b / 4)

    # --- fp32 Calculation (C) ---
    c = 0.0
    for num in sequence:
        c = c + num
    c = c * 16
    c = c + 0.25
    c = c / 4

    print(f"A (nf4) = {a}")
    print(f"B (bf16) = {b}")
    print(f"C (fp32) = {c}")
    print("\nThe final equation is: ceil((B - C - A) * 10)")
    print(f"Substituting the values: ceil(({b} - {c} - {a}) * 10)")

solve()