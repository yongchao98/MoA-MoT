import math

# Machine operation times in nanoseconds (ns)
int_add_time = 1
int_mul_time = 2
float_add_time = 9
float_mul_time = 19

# Size of the series
n = 1000

# --- Direct Convolution Calculation ---
print("Step 1: Calculating cost of Direct Convolution")
num_ops_direct = n**2
print(f"For n={n}, direct convolution requires approximately {num_ops_direct} additions and {num_ops_direct} multiplications.")
print("-" * 30)

# Direct convolution with integers
print("Calculation for Direct Convolution with Integers:")
time_direct_int_add = num_ops_direct * int_add_time
time_direct_int_mul = num_ops_direct * int_mul_time
total_time_direct_int = time_direct_int_add + time_direct_int_mul
print(f"Integer Additions Time: {num_ops_direct} * {int_add_time} ns = {time_direct_int_add} ns")
print(f"Integer Multiplications Time: {num_ops_direct} * {int_mul_time} ns = {time_direct_int_mul} ns")
print(f"Total Time = {time_direct_int_add} + {time_direct_int_mul} = {total_time_direct_int} ns")
print("-" * 30)

# Direct convolution with floating points
print("Calculation for Direct Convolution with Floating Points:")
time_direct_float_add = num_ops_direct * float_add_time
time_direct_float_mul = num_ops_direct * float_mul_time
total_time_direct_float = time_direct_float_add + time_direct_float_mul
print(f"Floating Point Additions Time: {num_ops_direct} * {float_add_time} ns = {time_direct_float_add} ns")
print(f"Floating Point Multiplications Time: {num_ops_direct} * {float_mul_time} ns = {time_direct_float_mul} ns")
print(f"Total Time = {time_direct_float_add} + {time_direct_float_mul} = {total_time_direct_float} ns")
print("-" * 30)

# --- FFT-based Convolution Calculation ---
print("\nStep 2: Calculating cost of FFT-based Convolution")
# Find the next power of 2 for FFT size
N = 2**(math.ceil(math.log2(2 * n - 1)))
log2_N = int(math.log2(N))

# Number of real operations for one FFT
# A complex multiplication = 4 real mult + 2 real add
# A complex addition = 2 real add
complex_mult_per_fft = (N / 2) * log2_N
complex_add_per_fft = N * log2_N
real_mult_per_fft = complex_mult_per_fft * 4
real_add_per_fft = complex_mult_per_fft * 2 + complex_add_per_fft * 2

# Total operations for 2 FFTs, 1 IFFT, and N element-wise complex multiplications
total_real_mult_fft = 3 * real_mult_per_fft + N * 4
total_real_add_fft = 3 * real_add_per_fft + N * 2

print(f"FFT size (N) must be a power of 2 >= {2*n-1}. Using N = {N}.")
print(f"Total real multiplications for FFT method: 3 * (({N}/2)*{log2_N} * 4) + ({N} * 4) = {int(total_real_mult_fft)}")
print(f"Total real additions for FFT method: 3 * (({N}/2)*{log2_N} * 2 + {N}*{log2_N} * 2) + ({N} * 2) = {int(total_real_add_fft)}")

# Calculate total time for FFT-based convolution
time_fft_mul = total_real_mult_fft * float_mul_time
time_fft_add = total_real_add_fft * float_add_time
total_time_fft = time_fft_mul + time_fft_add
print("\nCalculation for FFT-based Convolution Time:")
print(f"Floating Point Multiplications Time: {int(total_real_mult_fft)} * {float_mul_time} ns = {int(time_fft_mul)} ns")
print(f"Floating Point Additions Time: {int(total_real_add_fft)} * {float_add_time} ns = {int(time_fft_add)} ns")
print(f"Total Time = {int(time_fft_mul)} + {int(time_fft_add)} = {int(total_time_fft)} ns")
print("-" * 30)

# --- Comparison and Conclusion ---
print("\nStep 3: Comparing the total times")
print(f"Time for Direct Convolution (Integer): {total_time_direct_int} ns")
print(f"Time for FFT-based Convolution (Float): {int(total_time_fft)} ns")
print(f"Time for Direct Convolution (Float): {total_time_direct_float} ns")

times = {
    "Direct convolution with integers": total_time_direct_int,
    "FFT-based convolution": total_time_fft,
    "Direct convolution with floating points": total_time_direct_float
}

fastest_method = min(times, key=times.get)
print(f"\nConclusion: The fastest algorithm is '{fastest_method}'.")
