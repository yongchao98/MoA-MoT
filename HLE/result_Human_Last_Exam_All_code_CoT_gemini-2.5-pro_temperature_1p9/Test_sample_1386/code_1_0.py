import math
import struct

def to_bf16(num: float) -> float:
    """
    Simulates the conversion of a float to a bfloat16 by truncating the mantissa.
    A bf16 has the same exponent as fp32 but only 7 mantissa bits.
    This is equivalent to keeping the top 16 bits of a 32-bit float.
    """
    # Pack the float into its 32-bit IEE 754 representation
    f32_bytes = struct.pack('>f', num)
    # Convert bytes to an integer for bitwise operations
    int_val = int.from_bytes(f32_bytes, 'big')
    # Truncate by shifting right by 16 bits (removes the lower 16 bits),
    # then shifting back left to maintain the magnitude.
    bf16_int_val = (int_val >> 16) << 16
    # Convert the truncated integer back to a 4-byte sequence and unpack as a float
    return struct.unpack('>f', bf16_int_val.to_bytes(4, 'big'))[0]

def to_nf4(num: float) -> float:
    """
    Simulates the clamping behavior of the nf4 format to its specified range.
    """
    NF4_MIN = -8.0
    NF4_MAX = 7.5
    return max(NF4_MIN, min(num, NF4_MAX))

def run_simulation():
    """
    Runs the calculation sequence for fp32, bf16, and nf4 formats.
    """
    operations = [
        ('+', 7), ('+', 7), ('+', 0.125), ('-', 7), ('-', 7), 
        ('+', 7), ('+', 7), ('+', 0.0625)
    ]
    
    post_processing = [
        ('*', 16), ('+', 0.25), ('/', 4)
    ]

    # --- FP32 Simulation (Value C) ---
    c = 0.0
    for op, val in operations:
        c += val
    for op, val in post_processing:
        if op == '*': c *= val
        elif op == '+': c += val
        elif op == '/': c /= val

    # --- BF16 Simulation (Value B) ---
    b = 0.0
    for op, val in operations:
        b = to_bf16(b + val)
    for op, val in post_processing:
        if op == '*': b = to_bf16(b * val)
        elif op == '+': b = to_bf16(b + val)
        elif op == '/': b = to_bf16(b / val)
    
    # --- NF4 Simulation (Value A) ---
    a = 0.0
    for op, val in operations:
        a = to_nf4(a + val)
    for op, val in post_processing:
        if op == '*': a = to_nf4(a * val)
        elif op == '+': a = to_nf4(a + val)
        elif op == '/': a = to_nf4(a / val)

    # --- Print Results ---
    print("--- GPU Format Simulation Results ---")
    print(f"Final value for nf4 (A): {a}")
    print(f"Final value for bf16 (B): {b}")
    print(f"Final value for fp32 (C): {c}")
    print("-" * 35)

    # --- Final Calculation ---
    b_minus_c = b - c
    b_minus_c_minus_a = b_minus_c - a
    multiplied_by_10 = b_minus_c_minus_a * 10
    final_result = math.ceil(multiplied_by_10)

    print("\n--- Final Answer Calculation ---")
    print(f"Expression: ceil((B - C - A) * 10)")
    print(f"Step 1: Substitute values => ceil(({b} - {c} - {a}) * 10)")
    print(f"Step 2: B - C => ceil(({b_minus_c} - {a}) * 10)")
    print(f"Step 3: (B-C) - A => ceil(({b_minus_c_minus_a}) * 10)")
    print(f"Step 4: Multiply by 10 => ceil({multiplied_by_10})")
    print(f"Step 5: Ceiling => {final_result}")
    
    return final_result

# Run the full simulation and calculation
final_answer = run_simulation()
print(f"\n<<< {final_answer} >>>")
