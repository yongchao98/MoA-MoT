import struct
import math

def simulate_nf4():
    """
    Simulates the calculation using nf4 format.
    Assumes quantization to the nearest integer within a range of [-8, 7]
    to model the behavior of a 4-bit float with limited values.
    """
    # Define the quantization function for nf4
    def quantize_nf4(val):
        # Round to nearest integer
        rounded = round(val)
        # Clamp to the approximate integer range of nf4
        return max(-8.0, min(7.0, rounded))

    result = 0.0
    operations = [7, 7, 0.125, -7, -7, 7, 7, 0.0625]
    
    # Initial additions
    for op in operations:
        result += op
        result = quantize_nf4(result)

    # Final operations
    result *= 16
    result = quantize_nf4(result)
    
    result += 0.25
    result = quantize_nf4(result)
    
    result /= 4
    result = quantize_nf4(result)
    
    return result

def simulate_bf16():
    """
    Simulates the calculation using bf16 format.
    bf16 is simulated by truncating the mantissa of an fp32 number.
    """
    # Define the conversion function from a Python float (fp64) to a simulated bf16
    def to_bf16(f_val):
        # Pack float to its 32-bit integer representation
        try:
            int_val = struct.unpack('<I', struct.pack('<f', f_val))[0]
        except OverflowError:
            # Handle potential infinity
            return float('inf') if f_val > 0 else float('-inf')
            
        # Truncate the last 16 bits (the lower part of the mantissa)
        int_val &= 0xFFFF0000
        # Unpack back to a float
        return struct.unpack('<f', struct.pack('<I', int_val))[0]

    result = 0.0
    operations = [7, 7, 0.125, -7, -7, 7, 7, 0.0625]

    # Initial additions
    for op in operations:
        result += op
        result = to_bf16(result)

    # Final operations
    result *= 16
    result = to_bf16(result)
    
    result += 0.25
    result = to_bf16(result)
    
    result /= 4
    result = to_bf16(result)
    
    return result

def simulate_fp32():
    """
    Simulates the calculation using fp32 format.
    Standard Python floats (fp64) have enough precision for this.
    """
    result = 0.0
    operations = [7, 7, 0.125, -7, -7, 7, 7, 0.0625]

    # Initial additions
    for op in operations:
        result += op

    # Final operations
    result *= 16
    result += 0.25
    result /= 4
    
    return result

if __name__ == "__main__":
    # Calculate the values for A, B, and C
    A = simulate_nf4()
    B = simulate_bf16()
    C = simulate_fp32()

    # The problem asks to output each number in the final equation.
    # The final equation is ceil((B - C - A) * 10)
    # So we print the values for A, B, and C.
    
    print(f"A (nf4 result) = {A}")
    print(f"B (bf16 result) = {B}")
    print(f"C (fp32 result) = {C}")
    
    # The user is asked to perform the final calculation mentally.
    # For verification:
    # B - C - A = 56.75 - 56.8125 - 2.0 = -2.0625
    # (B - C - A) * 10 = -20.625
    # ceil(-20.625) = -20
    
    # The final answer is provided below as requested.