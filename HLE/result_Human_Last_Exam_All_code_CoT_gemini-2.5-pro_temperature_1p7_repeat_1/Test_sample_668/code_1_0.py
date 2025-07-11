import math

# --- Machine Operation Times in ns ---
int_add_time = 1
int_mul_time = 2
float_add_time = 9
float_mul_time = 19

# --- Problem Parameters ---
n = 1000

print("--- Calculation for Convolution Algorithms ---")

# --- 1. Direct Convolution ---
print("\n1. Direct Convolution Method (O(n^2))")
num_mul_direct = n * n
num_add_direct = n * n

# 1a. With Integers
time_direct_int = (num_mul_direct * int_mul_time) + (num_add_direct * int_add_time)
print("   - Integer Calculation:")
print(f"     Number of Multiplications = {num_mul_direct}")
print(f"     Number of Additions = {num_add_direct}")
print(f"     Estimated Time = ({num_mul_direct} * {int_mul_time} ns) + ({num_add_direct} * {int_add_time} ns) = {time_direct_int:,} ns")

# 1b. With Floating Points
time_direct_float = (num_mul_direct * float_mul_time) + (num_add_direct * float_add_time)
print("   - Floating Point Calculation:")
print(f"     Number of Multiplications = {num_mul_direct}")
print(f"     Number of Additions = {num_add_direct}")
print(f"     Estimated Time = ({num_mul_direct} * {float_mul_time} ns) + ({num_add_direct} * {float_add_time} ns) = {time_direct_float:,} ns")


# --- 2. FFT-based Convolution ---
print("\n2. FFT-based Convolution Method (O(N log N))")
# Padded size N must be a power of 2 and >= 2n-1
N = 1
while N < (2 * n - 1):
    N *= 2
log2N = int(math.log2(N))

# Operations for 1 complex FFT (Cooley-Tukey)
# A complex mult = 4 real mult + 2 real add
# A complex add = 2 real add
complex_mul_per_fft = (N / 2) * log2N
complex_add_per_fft = N * log2N
real_mul_per_fft = complex_mul_per_fft * 4
real_add_per_fft = (complex_mul_per_fft * 2) + (complex_add_per_fft * 2)

# Total operations for 2 forward FFTs, 1 inverse FFT, and 1 element-wise product
total_real_mul = (3 * real_mul_per_fft) + (N * 4) # 3 FFTs + N element-wise complex mult
total_real_add = (3 * real_add_per_fft) + (N * 2) # 3 FFTs + N element-wise complex mult

time_fft_float = (total_real_mul * float_mul_time) + (total_real_add * float_add_time)
print("   - Floating Point Calculation:")
print(f"     Padded size N = {N}")
print(f"     Total Floating Point Multiplications = {int(total_real_mul)}")
print(f"     Total Floating Point Additions = {int(total_real_add)}")
print(f"     Estimated Time = ({int(total_real_mul)} * {float_mul_time} ns) + ({int(total_real_add)} * {float_add_time} ns) = {int(time_fft_float):,} ns")

print("\n--- Conclusion ---")
if time_direct_int < time_fft_float and time_direct_int < time_direct_float:
    print("The fastest algorithm is Direct Convolution with Integers.")
elif time_fft_float < time_direct_int and time_fft_float < time_direct_float:
    print("The fastest algorithm is FFT-based Convolution with Floating Points.")
else:
    print("The fastest algorithm is Direct Convolution with Floating Points.")