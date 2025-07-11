import math

# --- Machine Operation Times ---
time_int_add = 1  # ns
time_int_mul = 2  # ns
time_float_add = 9 # ns
time_float_mul = 19 # ns

# --- Problem Parameters ---
n = 1000

print("--- Analysis of Convolution Algorithms ---")

# --- 1. Direct Convolution ---
print("\n1. Direct Convolution Calculation (Complexity O(n^2))")
ops_mul_direct = n * n
ops_add_direct = n * n
print(f"For n = {n}, the number of operations is roughly:")
print(f"- Multiplications: {n} * {n} = {ops_mul_direct}")
print(f"- Additions: {n} * {n} = {ops_add_direct}")

# B: Direct convolution with integers
time_direct_int_mul = ops_mul_direct * time_int_mul
time_direct_int_add = ops_add_direct * time_int_add
total_time_direct_int = time_direct_int_mul + time_direct_int_add
print("\nOption B: Direct Convolution with Integers")
print(f"Time = ({ops_mul_direct} multiplications * {time_int_mul} ns) + ({ops_add_direct} additions * {time_int_add} ns)")
print(f"Time = {time_direct_int_mul} ns + {time_direct_int_add} ns = {total_time_direct_int:,} ns")

# C: Direct convolution with floating points
time_direct_float_mul = ops_mul_direct * time_float_mul
time_direct_float_add = ops_add_direct * time_float_add
total_time_direct_float = time_direct_float_mul + time_direct_float_add
print("\nOption C: Direct Convolution with Floating Points")
print(f"Time = ({ops_mul_direct} multiplications * {time_float_mul} ns) + ({ops_add_direct} additions * {time_float_add} ns)")
print(f"Time = {time_direct_float_mul:,} ns + {time_direct_float_add:,} ns = {total_time_direct_float:,} ns")


# --- 2. FFT-based Convolution ---
print("\n2. FFT-based Convolution Calculation (Complexity O(N log N))")
output_len = 2 * n - 1
N = 1
while N < output_len:
    N *= 2
log2N = int(math.log2(N))
print(f"Padded FFT size N = {N} (next power of 2 >= {output_len}), log2(N) = {log2N}")

# Number of transforms: 2 forward FFTs, 1 inverse IFFT
num_transforms = 3

# Operations for 3 transforms
ops_mul_transforms = num_transforms * (2 * N * log2N)
ops_add_transforms = num_transforms * (3 * N * log2N)

# Operations for element-wise complex multiplication
ops_mul_pointwise = 4 * N
ops_add_pointwise = 2 * N

# Total operations for FFT method
total_ops_mul_fft = ops_mul_transforms + ops_mul_pointwise
total_ops_add_fft = ops_add_transforms + ops_add_pointwise

print(f"\nTotal operations for {num_transforms} transforms (FFT/IFFT):")
print(f"- Multiplications: {num_transforms} * (2 * {N} * {log2N}) = {ops_mul_transforms:,}")
print(f"- Additions: {num_transforms} * (3 * {N} * {log2N}) = {ops_add_transforms:,}")

print(f"\nTotal operations for pointwise complex product:")
print(f"- Multiplications: 4 * {N} = {ops_mul_pointwise:,}")
print(f"- Additions: 2 * {N} = {ops_add_pointwise:,}")

# A: FFT convolution
time_fft_mul = total_ops_mul_fft * time_float_mul
time_fft_add = total_ops_add_fft * time_float_add
total_time_fft = time_fft_mul + time_fft_add
print("\nOption A: FFT-based Convolution (Floating Point)")
print(f"Total multiplications = {ops_mul_transforms:,} + {ops_mul_pointwise:,} = {total_ops_mul_fft:,}")
print(f"Total additions = {ops_add_transforms:,} + {ops_add_pointwise:,} = {total_ops_add_fft:,}")
print(f"Time = ({total_ops_mul_fft:,} multiplications * {time_float_mul} ns) + ({total_ops_add_fft:,} additions * {time_float_add} ns)")
print(f"Time = {time_fft_mul:,} ns + {time_fft_add:,} ns = {total_time_fft:,} ns")

# --- 3. Conclusion ---
print("\n--- Conclusion ---")
print(f"Time for Direct Convolution (integer):      {total_time_direct_int:,} ns")
print(f"Time for FFT-based Convolution (float):   {total_time_fft:,} ns")
print(f"Time for Direct Convolution (float):      {total_time_direct_float:,} ns")

if total_time_direct_int < total_time_fft and total_time_direct_int < total_time_direct_float:
    fastest_algo = "B. Direct convolution with integers"
elif total_time_fft < total_time_direct_int and total_time_fft < total_time_direct_float:
    fastest_algo = "A. FFT"
else:
    fastest_algo = "C. Direct convolution with floating points"

print(f"\nThe fastest algorithm is: {fastest_algo}")