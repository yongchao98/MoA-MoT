import numpy as np

def quantize_nf4(val):
    """
    Simulates nf4 quantization.
    Based on the problem description, nf4 is modeled as a 4-bit format
    representing the 16 integers from -8 to 7.
    The simulation rounds the value to the nearest integer and then clips it
    to the representable range.
    """
    return np.clip(np.round(val), -8, 7)

def quantize_bf16(val):
    """
    Simulates bf16 rounding.
    bfloat16 has a 7-bit mantissa. The precision (step size between
    representable numbers) depends on the number's magnitude (exponent).
    """
    if val == 0:
        return 0.0
    # Calculate the exponent in base 2
    exponent = np.floor(np.log2(np.abs(val)))
    # The step size is 2^(exponent - mantissa_bits)
    step = 2**(exponent - 7)
    # Quantize the value by rounding to the nearest representable step
    return step * np.round(val / step)

def run_simulation(quantizer=None):
    """
    Runs the full calculation sequence, applying a specified quantizer
    after each arithmetic operation.
    """
    # The sequence of numbers to add
    sequence = [7, 7, 0.125, -7, -7, 7, 7, 0.0625]

    # Start with 0
    x = 0.0

    # Add numbers in sequence, quantizing after each step
    for num in sequence:
        x = x + num
        if quantizer:
            x = quantizer(x)

    # Perform final operations, quantizing after each step
    x = x * 16
    if quantizer:
        x = quantizer(x)

    x = x + 0.25
    if quantizer:
        x = quantizer(x)

    x = x / 4
    if quantizer:
        x = quantizer(x)

    return x

# --- Run Simulations ---

# A: nf4 simulation
val_a = run_simulation(quantizer=quantize_nf4)

# B: bf16 simulation
val_b = run_simulation(quantizer=quantize_bf16)

# C: fp32 simulation (no precision loss for these operations)
val_c = run_simulation(quantizer=None)

# --- Output Results ---
# The final equation requires the values A, B, and C.
print(f"Final value for nf4 (A) = {val_a}")
print(f"Final value for bf16 (B) = {val_b}")
print(f"Final value for fp32 (C) = {val_c}")