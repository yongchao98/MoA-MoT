import struct
import math

def quantize_nf4(value):
    """
    Simulates the nf4 format by rounding to the nearest integer and clamping.
    The range is [-8, 7.5] as specified.
    This captures the low precision and fixed range characteristics.
    """
    # Round to the nearest whole number to simulate low precision
    rounded_value = round(value)
    # Clamp the value to the specified range [-8, 7.5]
    clamped_value = max(-8.0, min(7.5, rounded_value))
    return clamped_value

def to_bf16(value):
    """
    Simulates a float32 to bfloat16 conversion.
    bfloat16 has the same 8-bit exponent as float32 but only a 7-bit mantissa.
    This is simulated by rounding the float32's 23-bit mantissa to 7 bits.
    """
    try:
        # Pack the float into its 32-bit integer representation
        packed_val = struct.pack('>f', value)
        int_val = struct.unpack('>I', packed_val)[0]
        
        # Add 0x8000 to handle rounding to nearest. This adds 1 to the
        # 16th bit, which is the first bit to be truncated. If it's 1,
        # the 17th bit (LSB of bf16 mantissa) will be incremented.
        int_val += 0x8000
        
        # Mask out the lower 16 bits of the mantissa to truncate it to 7 bits
        int_val &= 0xFFFF0000
        
        # Unpack the modified integer back into a float
        return struct.unpack('>f', struct.pack('>I', int_val))[0]
    except (OverflowError, ValueError):
        # Return original value if it's inf, nan, etc.
        return value

def run_simulation():
    """
    Runs the calculation sequence for all three number formats.
    """
    initial_ops = [
        ('add', 7), ('add', 7), ('add', 0.125), ('add', -7), 
        ('add', -7), ('add', 7), ('add', 7), ('add', 0.0625)
    ]
    final_ops = [
        ('mul', 16), ('add', 0.25), ('div', 4)
    ]

    # --- FP32 Calculation (C) ---
    # Standard Python floats have enough precision to be exact.
    c = 0.0
    for op, val in initial_ops:
        if op == 'add':
            c += val
    for op, val in final_ops:
        if op == 'mul':
            c *= val
        elif op == 'add':
            c += val
        elif op == 'div':
            c /= val

    # --- BF16 Calculation (B) ---
    # Each operation result is converted to bfloat16.
    b = 0.0
    for op, val in initial_ops:
        if op == 'add':
            b = to_bf16(b + val)
    for op, val in final_ops:
        if op == 'mul':
            b = to_bf16(b * val)
        elif op == 'add':
            b = to_bf16(b + val)
        elif op == 'div':
            b = to_bf16(b / val)

    # --- NF4 Calculation (A) ---
    # Each operation result is quantized to nf4.
    a = 0.0
    for op, val in initial_ops:
        if op == 'add':
            a = quantize_nf4(a + val)
    for op, val in final_ops:
        if op == 'mul':
            a = quantize_nf4(a * val)
        elif op == 'add':
            a = quantize_nf4(a + val)
        elif op == 'div':
            a = quantize_nf4(a / val)

    print("--- Simulation Results ---")
    print(f"The final value for nf4 (A) is: {a}")
    print(f"The final value for bf16 (B) is: {b}")
    print(f"The final value for fp32 (C) is: {c}")
    
    print("\n--- Final Equation ---")
    print("The equation to solve is: ceil((B - C - A) * 10)")
    print("Substituting the calculated values:")
    print(f"A = {a}")
    print(f"B = {b}")
    print(f"C = {c}")
    print(f"Equation: ceil(({b} - {c} - {a}) * 10)")

if __name__ == "__main__":
    run_simulation()