import numpy as np
import struct

def simulate_gpu_formats():
    """
    This function simulates a sequence of calculations on three different
    GPU number formats: nf4, bf16, and fp32.
    """
    sequence = [7, 7, 0.125, -7, -7, 7, 7, 0.0625]

    # --- Helper functions for quantization ---

    def quantize_nf4(value):
        """Simulates nf4 by clipping the value to its defined range."""
        return np.clip(value, -8.0, 7.5)

    def quantize_bf16(value):
        """
        Simulates bf16 rounding by taking a float32, representing it as 
        32 bits, masking the lower 16 bits of the mantissa, and converting back.
        """
        # Pack float to 4 bytes (32 bits)
        packed_float = struct.pack('!f', float(value))
        # Unpack as an integer to perform bitwise operations
        int_val = struct.unpack('!I', packed_float)[0]
        # Mask the lower 16 bits to truncate the mantissa
        int_val &= 0xFFFF0000
        # Pack back to bytes and unpack as a float
        unpacked_float = struct.unpack('!f', struct.pack('!I', int_val))[0]
        return unpacked_float

    def quantize_fp32(value):
        """Simulates fp32 by casting the value."""
        return np.float32(value)

    # --- Simulation loop ---
    
    formats = {
        'nf4': quantize_nf4,
        'bf16': quantize_bf16,
        'fp32': quantize_fp32
    }
    
    results = {}

    for name, quantize_func in formats.items():
        # Start with 0 in the correct format
        total = quantize_func(0)

        # Add numbers in sequence
        for num in sequence:
            total = quantize_func(total + num)
        
        # Perform post-processing
        total = quantize_func(total * 16)
        total = quantize_func(total + 0.25)
        total = quantize_func(total / 4)
        
        results[name] = total

    # --- Print results ---
    A = results['nf4']
    B = results['bf16']
    C = results['fp32']

    print("After running the simulation, we get the following final values:")
    print(f"Final value for nf4 (A): {A}")
    print(f"Final value for bf16 (B): {B}")
    print(f"Final value for fp32 (C): {C}")
    print("\nThe problem asks for the result of ceil((B - C - A) * 10)")
    print("The final equation with the computed values is:")
    print(f"ceil(({B} - {C} - {A}) * 10)")

simulate_gpu_formats()
