import struct
import numpy as np

def simulate_gpu_formats():
    """
    Simulates calculations on three different GPU number formats (nf4, bf16, fp32)
    and prints the final results.
    """

    # --- Helper functions for number formats ---

    def quantize_nf4(value):
        """
        Simulates the nf4 format by rounding to the nearest integer
        and then clipping to the 16-value range [-8, 7].
        """
        # Round to the nearest integer to simulate precision loss for fractions
        rounded = round(value)
        # Clip to the asymmetric range of 16 integer values [-8, 7]
        clipped = max(-8, min(7, rounded))
        return float(clipped)

    def to_bf16(value):
        """
        Simulates bfloat16 by converting a float to 32-bit representation,
        truncating it to 16 bits (bf16 format), and converting back.
        """
        # Pack the float into its 32-bit (single precision) byte representation
        # Use numpy.float32 to ensure we start from a 32-bit float
        f32_val = np.float32(value)
        packed = struct.pack('>f', f32_val)
        # Unpack as an integer to perform bitwise operations
        int_val, = struct.unpack('>I', packed)
        # Truncate by shifting right 16 bits, dropping the last 16 bits of the mantissa
        int_val = int_val >> 16
        # Shift back to make it a 32-bit representation again (padded with zeros)
        int_val = int_val << 16
        # Repack and unpack as a float
        bf16_val, = struct.unpack('>f', struct.pack('>I', int_val))
        return bf16_val

    # --- Initial values and sequence of operations ---
    
    val_nf4 = 0.0
    val_bf16 = 0.0
    val_fp32 = np.float32(0.0)

    operations = [7, 7, 0.125, -7, -7, 7, 7, 0.0625]

    # --- Perform sequence of additions ---
    
    for op in operations:
        val_nf4 = quantize_nf4(val_nf4 + op)
        val_bf16 = to_bf16(val_bf16 + op)
        val_fp32 = np.float32(val_fp32 + op)
    
    # --- Perform final operations ---
    
    # 1. Multiply by 16
    val_nf4 = quantize_nf4(val_nf4 * 16)
    val_bf16 = to_bf16(val_bf16 * 16)
    val_fp32 = np.float32(val_fp32 * 16)

    # 2. Add 0.25
    val_nf4 = quantize_nf4(val_nf4 + 0.25)
    val_bf16 = to_bf16(val_bf16 + 0.25)
    val_fp32 = np.float32(val_fp32 + 0.25)
    
    # 3. Divide by 4
    val_nf4 = quantize_nf4(val_nf4 / 4)
    val_bf16 = to_bf16(val_bf16 / 4)
    val_fp32 = np.float32(val_fp32 / 4)

    # --- Print the final values ---
    
    print(f"Final value for A (nf4) = {val_nf4}")
    print(f"Final value for B (bf16) = {val_bf16}")
    print(f"Final value for C (fp32) = {val_fp32}")

if __name__ == '__main__':
    simulate_gpu_formats()