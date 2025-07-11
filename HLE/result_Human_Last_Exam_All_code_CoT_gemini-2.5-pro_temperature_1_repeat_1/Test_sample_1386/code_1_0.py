import numpy as np
import math

# Define the sequence of numbers for the calculation
numbers = [7, 7, 0.125, -7, -7, 7, 7, 0.0625]

# --- Part 1: NF4 Simulation (Value A) ---
# nf4 is defined by its range [-8, 7.5]. Values are clipped if they exceed this range.
nf4_min, nf4_max = -8.0, 7.5
current_a = 0.0
for num in numbers:
    current_a += num
    # Clip the result of each addition to the nf4 range
    current_a = max(nf4_min, min(nf4_max, current_a))

# The subsequent operations are performed on the final clipped result
# using standard precision, as is typical for debugging analysis.
A = (current_a * 16 + 0.25) / 4

# --- Part 2: BF16 Simulation (Value B) ---
# bf16 has a 7-bit mantissa, which can lead to rounding errors.
# We use numpy.bfloat16 to simulate this accurately.
# Note: Requires numpy version 1.20 or newer.
bf16_dtype = np.bfloat16
current_b = bf16_dtype(0)
for num in numbers:
    current_b += bf16_dtype(num)

# Perform final operations, ensuring each step respects bf16 precision.
# The precision loss occurs when adding 0.25 to 227.
res_b = current_b * bf16_dtype(16)
res_b = res_b + bf16_dtype(0.25)
res_b = res_b / bf16_dtype(4)
B = float(res_b)

# --- Part 3: FP32 Simulation (Value C) ---
# fp32 has a 23-bit mantissa and can handle these calculations without precision loss.
fp32_dtype = np.float32
current_c = fp32_dtype(0)
for num in numbers:
    current_c += fp32_dtype(num)

# Perform final operations in fp32 precision.
res_c = current_c * fp32_dtype(16)
res_c = res_c + fp32_dtype(0.25)
res_c = res_c / fp32_dtype(4)
C = float(res_c)

# --- Final Output ---
# The prompt asks to output each number in the final equation.
print("The final calculation is: ceil((B - C - A) * 10)")
print(f"Value for A (from nf4) = {A}")
print(f"Value for B (from bf16) = {B}")
print(f"Value for C (from fp32) = {C}")

# The prompt asks for the final answer in a specific format after you
# have seen the values of A, B, and C.
# The calculation is: ceil((56.75 - 56.8125 - 30.0625) * 10)
# = ceil((-0.0625 - 30.0625) * 10)
# = ceil(-30.125 * 10)
# = ceil(-301.25)
# = -301