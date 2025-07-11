import math

# Machine operation times in ns
INT_ADD_TIME = 1
INT_MUL_TIME = 2
FLOAT_ADD_TIME = 9
FLOAT_MUL_TIME = 19

# Series length
n = 1000

print("--- Calculation for n = 1000 elements ---")

# --- Direct Convolution ---
direct_muls = n * n
direct_adds = n * (n - 1)

# Case B: Direct convolution with integers
time_direct_int = (direct_muls * INT_MUL_TIME) + (direct_adds * INT_ADD_TIME)
print("\nB. Direct convolution with integers:")
print(f"   Operations: {direct_muls} multiplications and {direct_adds} additions.")
print(f"   Time = {direct_muls} * {INT_MUL_TIME} ns + {direct_adds} * {INT_ADD_TIME} ns = {time_direct_int:,} ns")

# Case C: Direct convolution with floating points
time_direct_float = (direct_muls * FLOAT_MUL_TIME) + (direct_adds * FLOAT_ADD_TIME)
print("\nC. Direct convolution with floating points:")
print(f"   Operations: {direct_muls} multiplications and {direct_adds} additions.")
print(f"   Time = {direct_muls} * {FLOAT_MUL_TIME} ns + {direct_adds} * {FLOAT_ADD_TIME} ns = {time_direct_float:,} ns")

# --- FFT-based Convolution ---
# Find the next power of 2 for FFT size N, where N >= 2n - 1
N = 1
while N < 2 * n - 1:
    N *= 2
log2_N = math.log2(N)

# A complex multiplication costs 4 real muls and 2 real adds.
# A complex addition costs 2 real adds.
# Total operations for 2 forward FFTs, 1 inverse FFT, and 1 element-wise product.
# This is equivalent to 3 FFTs and 1 element-wise product in terms of complex operations.
total_complex_muls = 3 * (N / 2 * log2_N) + N
total_complex_adds = 3 * (N * log2_N)

# Convert complex operations to real floating-point operations
fft_muls = total_complex_muls * 4
fft_adds = total_complex_muls * 2 + total_complex_adds * 2

# Case A: FFT-based convolution
time_fft = (fft_muls * FLOAT_MUL_TIME) + (fft_adds * FLOAT_ADD_TIME)
print("\nA. FFT-based convolution:")
print(f"   FFT size N = {N}")
print(f"   Real Operations: {int(fft_muls)} multiplications and {int(fft_adds)} additions.")
print(f"   Time = {int(fft_muls)} * {FLOAT_MUL_TIME} ns + {int(fft_adds)} * {FLOAT_ADD_TIME} ns = {int(time_fft):,} ns")


# --- Comparison ---
results = {
    "FFT-based convolution (floats)": time_fft,
    "Direct convolution with integers": time_direct_int,
    "Direct convolution with floating points": time_direct_float,
}

fastest_method = min(results, key=results.get)
fastest_time = results[fastest_method]

print("\n--- Conclusion ---")
print(f"The fastest method is '{fastest_method}' with an estimated time of {int(fastest_time):,} ns.")
