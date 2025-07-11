import math

# Machine operation times in nanoseconds
int_add_time = 1
int_mul_time = 2
fp_add_time = 9
fp_mul_time = 19

# Size of the series
n = 1000

# --- 1. Direct Convolution Calculation ---

# Number of operations for direct convolution
direct_mults = n * n
direct_adds = n * (n - 1)

# Time for direct convolution with integers
time_direct_int = (direct_mults * int_mul_time) + (direct_adds * int_add_time)

# Time for direct convolution with floating points
time_direct_fp = (direct_mults * fp_mul_time) + (direct_adds * fp_add_time)

print("--- Calculation Details ---")

# Print Integer Direct Convolution Calculation
print("\nDirect Convolution (Integers):")
print(f"({direct_mults} multiplications * {int_mul_time} ns) + ({direct_adds} additions * {int_add_time} ns) = {time_direct_int} ns")

# Print Floating Point Direct Convolution Calculation
print("\nDirect Convolution (Floating Points):")
print(f"({direct_mults} multiplications * {fp_mul_time} ns) + ({direct_adds} additions * {fp_add_time} ns) = {time_direct_fp} ns")

# --- 2. FFT-based Convolution Calculation ---

# Find the next power of 2 for FFT size (N >= 2n-1)
N = 1
while N < (2 * n - 1):
    N *= 2
log2_N = math.log2(N)

# Operations for one complex FFT
# Complex mult = 4 real mults + 2 real adds
# Complex add = 2 real adds
# One FFT has (N/2)log2(N) complex mults and N*log2(N) complex adds
fft_mults_one = 4 * (N / 2) * log2_N
# Note: simplifying 2*(N/2)log2N + 2*N*log2N = N*log2N + 2*N*log2N = 3*N*log2N
fft_adds_one = 3 * N * log2_N

# Total operations for FFT convolution:
# 2 FFTs, 1 IFFT (same cost), N complex point-wise multiplications
total_fft_mults = 3 * fft_mults_one + (N * 4) # 3 FFTs + N complex mults
total_fft_adds = 3 * fft_adds_one + (N * 2) # 3 FFTs + N complex mults

# Time for FFT-based convolution (always floating point)
time_fft = (total_fft_mults * fp_mul_time) + (total_fft_adds * fp_add_time)

# Print FFT Convolution Calculation
print("\nFFT-based Convolution (Floating Points):")
print(f"({int(total_fft_mults)} multiplications * {fp_mul_time} ns) + ({int(total_fft_adds)} additions * {fp_add_time} ns) = {int(time_fft)} ns")

# --- 3. Conclusion ---
print("\n--- Conclusion ---")
if time_direct_int < time_fft and time_direct_int < time_direct_fp:
    print("The fastest algorithm is Direct convolution with integers.")
elif time_fft < time_direct_int and time_fft < time_direct_fp:
    print("The fastest algorithm is FFT-based convolution.")
else:
    print("The fastest algorithm is Direct convolution with floating points.")
