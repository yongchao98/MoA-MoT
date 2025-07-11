import struct
import math

def get_gpu_format_results():
    """
    Calculates the result of a sequence of operations under three different
    simulated GPU number formats: nf4, bf16, and fp32.
    """

    # Helper function to simulate bfloat16 rounding.
    # It converts a 32-bit float to 16-bit bfloat by clearing the last 16 bits of the mantissa.
    def to_bfloat16(f_val):
        # Packs the float into 4 bytes (big-endian), truncates the last 2 bytes, and unpacks it.
        # This effectively simulates the precision loss of bfloat16.
        packed_fp32 = struct.pack('>f', f_val)
        truncated = packed_fp32[:2] + b'\x00\x00'
        return struct.unpack('>f', truncated)[0]

    # The sequence of numbers for the initial calculation
    numbers = [7, 7, 0.125, -7, -7, 7, 7, 0.0625]

    # --- Case 1: nf4 Simulation (A) ---
    # nf4 is simulated by clamping the result of each addition to its defined range [-8, 7.5].
    val_nf4 = 0.0
    for num in numbers:
        val_nf4 += num
        # Clamp the value to the nf4 range
        if val_nf4 > 7.5:
            val_nf4 = 7.5
        elif val_nf4 < -8.0:
            val_nf4 = -8.0

    # Final operations on the nf4 result
    # These are assumed to be done in higher precision, so no extra clamping here.
    A = ((val_nf4 * 16) + 0.25) / 4

    # --- Case 2: bfloat16 Simulation (B) ---
    # The intermediate sum is calculated first. For these numbers, no precision
    # is lost in the initial sum as all are perfectly representable.
    # val_bf16 = sum(numbers) -> this is 14.1875 which is exact in bf16
    val_bf16 = 14.1875

    # The final operations are performed with bfloat16 rounding at each step.
    res_b = to_bfloat16(val_bf16 * 16)  # 227.0, which is exact
    res_b = to_bfloat16(res_b + 0.25)   # 227.25 rounds to 227.0 in bf16
    B = to_bfloat16(res_b / 4)          # 227.0 / 4 = 56.75, which is exact in bf16

    # --- Case 3: fp32 Simulation (C) ---
    # Standard Python floats (double precision) are used, which is more than
    # enough to represent fp32 without error for this problem.
    val_fp32 = sum(numbers)

    # Final operations on the fp32 result
    C = ((val_fp32 * 16) + 0.25) / 4

    # --- Output Results ---
    print("Final values calculated:")
    print(f"A (nf4) = {A}")
    print(f"B (bf16) = {B}")
    print(f"C (fp32) = {C}")
    print("\nFinal Equation:")
    print(f"B - C - A = {B} - {C} - {A}")
    
    # Mental calculation:
    # B - C - A = 56.75 - 56.8125 - 30.0625
    #           = -0.0625 - 30.0625
    #           = -30.125
    # (B - C - A) * 10 = -301.25
    # ceil(-301.25) = -301

if __name__ == '__main__':
    get_gpu_format_results()