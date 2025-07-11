import math
import struct

def calculate_nf4():
    """
    Simulates the calculation in nf4 format.
    The key feature of nf4 is its limited range, causing clamping.
    Range: -8 to 7.5
    """
    print("--- Simulating nf4 (Value A) ---")
    
    val = 0.0
    min_val, max_val = -8.0, 7.5
    numbers = [7, 7, 0.125, -7, -7, 7, 7, 0.0625]
    
    # Perform summation with clamping after each step
    for num in numbers:
        val += num
        # Clamp the value to the nf4 range
        val = max(min_val, min(max_val, val))

    nf4_sum = val
    print(f"Sum after nf4 clamping: {nf4_sum}")
    
    # The final operations are performed on the result
    A = (nf4_sum * 16 + 0.25) / 4
    
    print(f"Final Equation for A: ({nf4_sum} * 16 + 0.25) / 4")
    print(f"A = {A}\n")
    return A

def quantize_bf16(n):
    """
    Quantizes a float to bfloat16 precision.
    This is done by taking the 32-bit representation of the float
    and zeroing out the last 16 bits, which truncates the mantissa.
    """
    # Pack the float into its 32-bit single-precision representation
    s = struct.pack('>f', n)
    # Unpack as a 32-bit integer
    i = struct.unpack('>I', s)[0]
    # Zero out the lower 16 bits (the end of the fp32 mantissa)
    i &= 0xFFFF0000
    # Pack the modified integer back into bytes
    s = struct.pack('>I', i)
    # Unpack back to a float
    return struct.unpack('>f', s)[0]

def calculate_bf16():
    """
    Simulates the calculation in bf16 format.
    The key feature is precision loss after every operation.
    """
    print("--- Simulating bf16 (Value B) ---")
    
    val = 0.0
    numbers = [7, 7, 0.125, -7, -7, 7, 7, 0.0625]

    # Perform summation with quantization after each step
    for num in numbers:
        val = quantize_bf16(val + num)
        
    bf16_sum = val
    print(f"Sum after bf16 quantization: {bf16_sum}")
    
    # Perform final operations with quantization after each step
    print(f"Step 1: {bf16_sum} * 16")
    val = quantize_bf16(val * 16)
    print(f"  Result after quantization: {val}")
    
    print(f"Step 2: {val} + 0.25")
    val = quantize_bf16(val + 0.25)
    print(f"  Result after quantization: {val}")

    print(f"Step 3: {val} / 4")
    val = quantize_bf16(val / 4)
    print(f"  Result after quantization: {val}")
    
    B = val
    print(f"B = {B}\n")
    return B

def calculate_fp32():
    """
    Simulates the calculation in fp32 format.
    Python's standard floats have enough precision for an exact simulation.
    """
    print("--- Simulating fp32 (Value C) ---")

    val = 0.0
    numbers = [7, 7, 0.125, -7, -7, 7, 7, 0.0625]

    for num in numbers:
        val += num
        
    fp32_sum = val
    print(f"Sum in fp32 precision: {fp32_sum}")
    
    C = (fp32_sum * 16 + 0.25) / 4
    
    print(f"Final Equation for C: ({fp32_sum} * 16 + 0.25) / 4")
    print(f"C = {C}\n")
    return C

if __name__ == "__main__":
    A = calculate_nf4()
    B = calculate_bf16()
    C = calculate_fp32()
    
    print("-----------------------------------------")
    print("Final Values:")
    print(f"A (nf4)  = {A}")
    print(f"B (bf16) = {B}")
    print(f"C (fp32) = {C}")
    print("-----------------------------------------")
    print("Now, mentally calculate ceil((B - C - A) * 10) to find the answer.")
