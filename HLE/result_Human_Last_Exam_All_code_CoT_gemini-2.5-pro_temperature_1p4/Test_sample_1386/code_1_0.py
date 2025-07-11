import struct
import math

def simulate_gpu_formats():
    """
    Simulates a sequence of floating-point operations on three different
    GPU number formats: nf4, bf16, and fp32.
    """
    
    # The sequence of numbers to add
    numbers = [7, 7, 0.125, -7, -7, 7, 7, 0.0625]

    # --- 1. nf4 Simulation (Value A) ---
    # nf4 has a limited range [-8, 7.5]. We simulate by clamping.
    x_nf4 = 0.0
    nf4_min = -8.0
    nf4_max = 7.5

    for num in numbers:
        x_nf4 += num
        # Clamp the result to the nf4 range
        if x_nf4 > nf4_max:
            x_nf4 = nf4_max
        elif x_nf4 < nf4_min:
            x_nf4 = nf4_min
    
    # Apply final operations for nf4
    A = (x_nf4 * 16 + 0.25) / 4

    # --- 2. bf16 Simulation (Value B) ---
    # bf16 has a 7-bit mantissa. We simulate by quantizing after each step.
    
    def quantize_to_bf16(f):
        """Quantizes an fp32 float to bf16 by truncating the mantissa."""
        try:
            # Pack float to its 32-bit integer representation
            s = struct.pack('>f', f)
            i = struct.unpack('>I', s)[0]
            # Shift right to discard the lower 16 bits (bf16 has 16 bits total)
            i = i >> 16
            # Shift back left to form a 32-bit number with zeroed-out lower 16 bits
            i = i << 16
            # Unpack back to a float
            s = struct.pack('>I', i)
            return struct.unpack('>f', s)[0]
        except (OverflowError, struct.error):
            return f # Return original if it's inf/nan

    x_bf16 = 0.0
    for num in numbers:
        x_bf16 += num
        x_bf16 = quantize_to_bf16(x_bf16)
    
    # Apply final operations for bf16, quantizing at each step
    x_bf16 = quantize_to_bf16(x_bf16 * 16)
    x_bf16 = quantize_to_bf16(x_bf16 + 0.25)
    B = quantize_to_bf16(x_bf16 / 4)

    # --- 3. fp32 Simulation (Value C) ---
    # fp32 is the standard Python float. No special handling is needed.
    x_fp32 = 0.0
    for num in numbers:
        x_fp32 += num
        
    # Apply final operations for fp32
    C = (x_fp32 * 16 + 0.25) / 4
    
    # --- Print Results ---
    print("--- GPU Number Format Detection Challenge ---")
    print(f"Sequence of additions: {numbers}")
    print("Final operations: (result * 16 + 0.25) / 4")
    print("-" * 45)
    print(f"Value A (nf4): {A}")
    print(f"Value B (bf16): {B}")
    print(f"Value C (fp32): {C}")
    print("-" * 45)
    print("Final equation to solve:")
    # The prompt asks to output each number in the final equation
    print(f"ceil(({B} - {C} - {A}) * 10)")

simulate_gpu_formats()