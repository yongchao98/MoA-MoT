import math
import struct

def simulate_gpu_formats():
    """
    Simulates a sequence of calculations on three different GPU number formats:
    nf4, bf16, and fp32 to find their final resulting values A, B, and C.
    """

    # --- Part 1: NF4 Simulation (Value A) ---
    def quantize_nf4(val):
        """Simulates nf4 by clamping the value to its specified range."""
        return max(-8.0, min(val, 7.5))

    val_a = 0.0
    sequence = [7, 7, 0.125, -7, -7, 7, 7, 0.0625]
    for num in sequence:
        val_a = quantize_nf4(val_a + num)
    
    val_a = quantize_nf4(val_a * 16)
    val_a = quantize_nf4(val_a + 0.25)
    val_a = quantize_nf4(val_a / 4)
    A = val_a

    # --- Part 2: BF16 Simulation (Value B) ---
    def quantize_bf16(f_val):
        """
        Simulates bfloat16 quantization by truncating the mantissa of a 32-bit float.
        This is done by taking the first 2 bytes of the fp32 representation.
        """
        # Pack the float into its 4-byte (big-endian) fp32 representation
        packed_fp32 = struct.pack('>f', f_val)
        # Create the bfloat16 representation by taking the first 2 bytes
        # and padding with two zero bytes.
        packed_bf16 = packed_fp32[:2] + b'\x00\x00'
        # Unpack back to a float
        return struct.unpack('>f', packed_bf16)[0]

    val_b = 0.0
    for num in sequence:
        # For bf16, intermediate results are exact for this specific sequence
        val_b += num

    val_b = quantize_bf16(val_b * 16)
    val_b = quantize_bf16(val_b + 0.25) # Precision loss occurs here
    val_b = quantize_bf16(val_b / 4)
    B = val_b
    
    # --- Part 3: FP32 Simulation (Value C) ---
    # Python's standard floats (fp64) are more than precise enough to model fp32
    val_c = 0.0
    for num in sequence:
        val_c += num
        
    val_c *= 16
    val_c += 0.25
    val_c /= 4
    C = val_c
    
    # --- Final Output ---
    print(f"Final Value for nf4 (A): {A}")
    print(f"Final Value for bf16 (B): {B}")
    print(f"Final Value for fp32 (C): {C}")
    print("\n--- Final Equation ---")
    print(f"The equation is: ceil(({B} - {C} - {A}) * 10)")

    # Perform the final calculation based on the computed values
    inner_val = (B - C - A) * 10
    final_result = math.ceil(inner_val)
    print(f"The result is: ceil({inner_val:.4f})")
    print(f"The final integer answer is: {final_result}")


simulate_gpu_formats()