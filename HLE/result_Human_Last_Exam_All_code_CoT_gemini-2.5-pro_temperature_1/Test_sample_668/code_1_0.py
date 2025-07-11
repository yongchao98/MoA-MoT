import math

# Machine operation times in nanoseconds (ns)
time_int_add = 1
time_int_mul = 2
time_float_add = 9
time_float_mul = 19

# Input series size
n = 1000

print("--- Analysis of Convolution Algorithms ---")
print(f"For two series of length n = {n}\n")

# --- Direct Convolution Analysis ---
print("Step 1: Analyzing Direct Convolution (O(n^2) complexity)")
# For simplicity, we approximate the number of operations as n^2
num_mult_direct = n * n
num_add_direct = n * n
print(f"This method requires approximately {num_mult_direct} multiplications and {num_add_direct} additions.")
print("-" * 40)

# Case B: Direct convolution with integers
print("Case B: Direct convolution with integers")
time_b_mult = num_mult_direct * time_int_mul
time_b_add = num_add_direct * time_int_add
total_time_b = time_b_mult + time_b_add
print(f"Time for integer multiplications = {num_mult_direct} * {time_int_mul} ns = {time_b_mult} ns")
print(f"Time for integer additions = {num_add_direct} * {time_int_add} ns = {time_b_add} ns")
print(f"Total estimated time = {time_b_mult} + {time_b_add} = {total_time_b} ns")
print("")

# Case C: Direct convolution with floating points
print("Case C: Direct convolution with floating points")
time_c_mult = num_mult_direct * time_float_mul
time_c_add = num_add_direct * time_float_add
total_time_c = time_c_mult + time_c_add
print(f"Time for floating point multiplications = {num_mult_direct} * {time_float_mul} ns = {time_c_mult} ns")
print(f"Time for floating point additions = {num_add_direct} * {time_float_add} ns = {time_c_add} ns")
print(f"Total estimated time = {time_c_mult} + {time_c_add} = {total_time_c} ns")
print("\n")


# --- FFT-based Convolution Analysis ---
print("Step 2: Analyzing FFT-based Convolution (O(N log N) complexity)")
# Determine FFT size N
output_len = 2 * n - 1
N = 1
while N < output_len:
    N *= 2
log2_N = int(math.log2(N))
print(f"FFT size N is the next power of 2 >= (2*n-1), so N = {N}.")
print("-" * 40)

# Calculate number of complex operations for 3 FFTs and 1 pointwise multiplication
num_complex_mult_3_ffts = 3 * (N / 2) * log2_N
num_complex_add_3_ffts = 3 * N * log2_N
num_complex_mult_pointwise = N
total_complex_mult = num_complex_mult_3_ffts + num_complex_mult_pointwise
total_complex_add = num_complex_add_3_ffts

# Convert to real floating-point operations
total_real_mult = total_complex_mult * 4
total_real_add = (total_complex_mult * 2) + (total_complex_add * 2)

print("Case A: FFT-based convolution with floating points")
print(f"Total real multiplications = {int(total_complex_mult)} complex mults * 4 = {int(total_real_mult)}")
print(f"Total real additions = ({int(total_complex_mult)} * 2) + ({int(total_complex_add)} * 2) = {int(total_real_add)}")

time_a_mult = total_real_mult * time_float_mul
time_a_add = total_real_add * time_float_add
total_time_a = time_a_mult + time_a_add
print(f"Time for floating point multiplications = {int(total_real_mult)} * {time_float_mul} ns = {int(time_a_mult)} ns")
print(f"Time for floating point additions = {int(total_real_add)} * {time_float_add} ns = {int(time_a_add)} ns")
print(f"Total estimated time = {int(time_a_mult)} + {int(time_a_add)} = {int(total_time_a)} ns")
print("\n")

# --- Conclusion ---
print("Step 3: Conclusion")
print("-" * 40)
print(f"Time for A (FFT) = {int(total_time_a)} ns")
print(f"Time for B (Direct Integer) = {total_time_b} ns")
print(f"Time for C (Direct Float) = {total_time_c} ns")
print("")

times = {
    "A": total_time_a,
    "B": total_time_b,
    "C": total_time_c
}

fastest_method_key = min(times, key=times.get)
fastest_method_name = ""
if fastest_method_key == 'A':
    fastest_method_name = "FFT"
elif fastest_method_key == 'B':
    fastest_method_name = "Direct convolution with integers"
elif fastest_method_key == 'C':
    fastest_method_name = "Direct convolution with floating points"

print(f"The fastest algorithm is B: {fastest_method_name}.")
