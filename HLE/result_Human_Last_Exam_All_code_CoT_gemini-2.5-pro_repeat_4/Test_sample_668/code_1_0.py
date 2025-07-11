import math

# Machine operation times in nanoseconds
int_add_time = 1
int_mul_time = 2
fp_add_time = 9
fp_mul_time = 19

# Series length
n = 1000

print("Step-by-step calculation to find the fastest convolution algorithm:")
print("-" * 60)

# --- 1. Direct Convolution ---
print("1. Direct Convolution Analysis (Complexity O(n^2))")
num_mults_direct = n**2
num_adds_direct = n**2
print(f"   For n={n}, this requires approximately {num_mults_direct:,} multiplications and {num_adds_direct:,} additions.")

# a) With Integers
time_direct_int_mult = num_mults_direct * int_mul_time
time_direct_int_add = num_adds_direct * int_add_time
total_time_direct_int = time_direct_int_mult + time_direct_int_add
print("\n   a) Using Integers:")
print(f"      Multiplication time = {n}^2 * {int_mul_time} ns = {time_direct_int_mult:,} ns")
print(f"      Addition time       = {n}^2 * {int_add_time} ns = {time_direct_int_add:,} ns")
print(f"      Total Integer Time  = {time_direct_int_mult:,} + {time_direct_int_add:,} = {total_time_direct_int:,} ns")

# b) With Floating Points
time_direct_fp_mult = num_mults_direct * fp_mul_time
time_direct_fp_add = num_adds_direct * fp_add_time
total_time_direct_fp = time_direct_fp_mult + time_direct_fp_add
print("\n   b) Using Floating Points:")
print(f"      Multiplication time = {n}^2 * {fp_mul_time} ns = {time_direct_fp_mult:,} ns")
print(f"      Addition time       = {n}^2 * {fp_add_time} ns = {time_direct_fp_add:,} ns")
print(f"      Total FP Time       = {time_direct_fp_mult:,} + {time_direct_fp_add:,} = {total_time_direct_fp:,} ns")
print("-" * 60)

# --- 2. FFT-based Convolution ---
print("2. FFT-based Convolution Analysis (Complexity O(N log N))")
# Find the next power of 2 for FFT size N >= 2n-1
N = 1
while N < 2 * n - 1:
    N *= 2
log2_N = int(math.log2(N))
print(f"   Padded size for FFT: N = {N} (since N >= 2*{n}-1 = {2*n-1})")
print(f"   log2({N}) = {log2_N}")

# Number of complex operations
# 2 FFTs + 1 IFFT = 3 transforms
# Point-wise multiplication of N elements
num_complex_mults_fft = 3 * (N / 2) * log2_N
num_complex_adds_fft = 3 * (N) * log2_N
num_complex_mults_pointwise = N
total_complex_mults = num_complex_mults_fft + num_complex_mults_pointwise
total_complex_adds = num_complex_adds_fft

# Convert to real FP operations
# 1 complex mult = 4 real mults + 2 real adds
# 1 complex add = 2 real adds
num_fp_mults_fft = total_complex_mults * 4
num_fp_adds_from_mult = total_complex_mults * 2
num_fp_adds_from_add = total_complex_adds * 2
total_fp_mults = int(num_fp_mults_fft)
total_fp_adds = int(num_fp_adds_from_mult + num_fp_adds_from_add)

print("\n   Calculating total real floating point operations:")
print(f"      Total complex multiplications = 3 * ({N}/2)*{log2_N} + {N} = {int(total_complex_mults):,}")
print(f"      Total complex additions       = 3 * {N}*{log2_N} = {int(total_complex_adds):,}")
print(f"      Total real multiplications = {int(total_complex_mults):,} * 4 = {total_fp_mults:,}")
print(f"      Total real additions       = ({int(total_complex_mults):,} * 2) + ({int(total_complex_adds):,} * 2) = {total_fp_adds:,}")

# Calculate total time for FFT method
time_fft_fp_mult = total_fp_mults * fp_mul_time
time_fft_fp_add = total_fp_adds * fp_add_time
total_time_fft = time_fft_fp_mult + time_fft_fp_add

print("\n   Calculating total time for FFT-based method:")
print(f"      Multiplication time = {total_fp_mults:,} * {fp_mul_time} ns = {time_fft_fp_mult:,} ns")
print(f"      Addition time       = {total_fp_adds:,} * {fp_add_time} ns = {time_fft_fp_add:,} ns")
print(f"      Total FFT Time      = {time_fft_fp_mult:,} + {time_fft_fp_add:,} = {total_time_fft:,} ns")
print("-" * 60)

# --- 3. Conclusion ---
print("3. Conclusion")
print(f"   Time for Direct Convolution (Integer):       {int(total_time_direct_int):,} ns")
print(f"   Time for FFT-based Convolution (FP):         {int(total_time_fft):,} ns")
print(f"   Time for Direct Convolution (Floating Point): {int(total_time_direct_fp):,} ns")

if total_time_direct_int < total_time_fft and total_time_direct_int < total_time_direct_fp:
    fastest_method = "Direct convolution with integers"
elif total_time_fft < total_time_direct_fp:
    fastest_method = "FFT-based convolution"
else:
    fastest_method = "Direct convolution with floating points"

print(f"\n   The fastest algorithm is: {fastest_method}.")