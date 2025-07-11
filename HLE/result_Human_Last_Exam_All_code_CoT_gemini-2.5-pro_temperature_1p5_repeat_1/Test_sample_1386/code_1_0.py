import numpy as np

def solve_gpu_challenge():
    """
    Simulates a sequence of floating-point operations on three different
    GPU number formats: nf4, bf16, and fp32.
    """
    
    # The sequence of numbers to add
    numbers_to_add = [7, 7, 0.125, -7, -7, 7, 7, 0.0625]

    # --- 1. nf4 Simulation (Value A) ---
    def simulate_nf4():
        # Define the 16 discrete values for nf4 based on the problem's
        # asymmetric range [-8, 7.5].
        nf4_values = np.array([-8., -7., -6., -5., -4., -3., -2., -1., 0., 1., 2., 3., 4., 5., 6., 7.5])
        
        def quantize_nf4(x):
            """Clips and rounds a value to the nearest nf4 value."""
            # Clip the value to the defined range of nf4
            clipped_x = np.clip(x, -8.0, 7.5)
            # Find the index of the closest value in our nf4 set
            closest_idx = np.argmin(np.abs(nf4_values - clipped_x))
            return nf4_values[closest_idx]

        acc = 0.0
        # print("nf4 trace:")
        # print(f"  start: {acc}")
        for num in numbers_to_add:
            acc += num
            acc = quantize_nf4(acc)
            # print(f"  + {num:<5} -> {acc}")

        # Final operations with quantization at each step
        acc *= 16
        acc = quantize_nf4(acc)
        # print(f"  * 16    -> {acc}")
        
        acc += 0.25
        acc = quantize_nf4(acc)
        # print(f"  + 0.25  -> {acc}")

        acc /= 4
        acc = quantize_nf4(acc)
        # print(f"  / 4     -> {acc}")
        
        return acc

    # --- 2. bf16 Simulation (Value B) ---
    def simulate_bf16():
        
        def quantize_bf16(f_val: np.float32) -> np.float32:
            """
            Simulates rounding a float32 value to bfloat16 precision.
            A bfloat16 has a 7-bit mantissa. This is achieved by taking a 
            32-bit float's integer representation and clearing the lower 16 bits.
            We add 0x8000 before clearing to implement round-to-nearest.
            """
            int_val = f_val.view(np.uint32)
            # Add rounding offset and then truncate the lower 16 bits
            int_val += 0x8000
            int_val &= 0xFFFF0000
            return int_val.view(np.float32)

        acc = np.float32(0.0)
        for num in numbers_to_add:
            acc = quantize_bf16(acc + np.float32(num))

        # Final operations with quantization at each step
        acc = quantize_bf16(acc * np.float32(16))
        acc = quantize_bf16(acc + np.float32(0.25))
        acc = quantize_bf16(acc / np.float32(4))
        
        return float(acc)

    # --- 3. fp32 Simulation (Value C) ---
    def simulate_fp32():
        acc = np.float32(0.0)
        for num in numbers_to_add:
            acc += np.float32(num)
            
        # Final operations (no quantization needed as precision is sufficient)
        acc *= np.float32(16)
        acc += np.float32(0.25)
        acc /= np.float32(4)

        return float(acc)

    # Calculate the final values for each format
    val_a = simulate_nf4()
    val_b = simulate_bf16()
    val_c = simulate_fp32()

    # The problem asks to output the equation with the final numbers
    print(f"The final calculated values are:")
    print(f"A (nf4)  = {val_a}")
    print(f"B (bf16) = {val_b}")
    print(f"C (fp32) = {val_c}")
    print("\nTo get the final answer, calculate ceil((B - C - A) * 10)")
    print(f"Based on the values above: ceil(({val_b} - {val_c} - {val_a}) * 10)")


if __name__ == '__main__':
    solve_gpu_challenge()
