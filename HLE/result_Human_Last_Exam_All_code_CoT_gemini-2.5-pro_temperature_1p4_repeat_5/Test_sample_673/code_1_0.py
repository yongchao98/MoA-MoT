import math

# Step 1: Define the given parameters
n = 1000  # Typical vector size
time_float_op = 5  # ns per floating point operation
time_int_op = 1    # ns per integer operation
time_call = 15     # ns per function call

print("--- Analysis of Algorithm 1: FFT-based ---")

# Step 2: Calculate the total time for the FFT-based algorithm (T_fft)

# For a linear convolution of two vectors of size n, the result is of size 2n-1.
# The FFT size N must be at least 2n-1. We choose the next power of 2.
required_len = 2 * n - 1
N = 2**math.ceil(math.log2(required_len))

# A convolution requires 3 FFTs (2 forward, 1 inverse).
# A recursive FFT of size N has (2N-2) function calls.
num_ffts = 3
num_calls_per_fft = 2 * N - 2
total_calls = num_ffts * num_calls_per_fft

# Calculate the time for function calls
time_for_calls = total_calls * time_call

# The problem states the algorithm performs 4n floating point operations.
num_float_ops_fft = 4 * n
time_for_float_ops_fft = num_float_ops_fft * time_float_op

# Calculate the total time for the FFT algorithm
T_fft = time_for_calls + time_for_float_ops_fft

print(f"Vector size n = {n}")
print(f"Required convolution length = 2*n - 1 = {required_len}")
print(f"FFT size N (next power of 2) = {N}\n")

print("Calculating FFT-based algorithm time:")
print(f"Total function calls = {num_ffts} FFTs * (2*{N} - 2) calls/FFT = {total_calls} calls")
print(f"Time for calls = {total_calls} * {time_call} ns/call = {time_for_calls} ns")
print(f"Floating point operations = 4 * n = {num_float_ops_fft} ops")
print(f"Time for float ops = {num_float_ops_fft} * {time_float_op} ns/op = {time_for_float_ops_fft} ns\n")
print(f"Total Time (FFT) = Time for calls + Time for float ops")
print(f"Total Time (FFT) = {time_for_calls} ns + {time_for_float_ops_fft} ns = {T_fft} ns\n")


print("--- Analysis of Algorithm 2: Direct Integer Convolution ---")

# Step 3: Calculate the total time for the direct integer algorithm (T_direct)

# The algorithm requires 2n floating point operations for conversion.
num_float_ops_direct = 2 * n
time_for_conversion = num_float_ops_direct * time_float_op

# The algorithm requires 2n^2 integer operations for convolution.
num_int_ops_direct = 2 * n**2
time_for_convolution = num_int_ops_direct * time_int_op

# Calculate the total time for the direct algorithm
T_direct = time_for_conversion + time_for_convolution

print("Calculating direct integer algorithm time:")
print(f"Floating point operations for conversion = 2 * n = {num_float_ops_direct} ops")
print(f"Time for conversion = {num_float_ops_direct} * {time_float_op} ns/op = {time_for_conversion} ns")
print(f"Integer operations for convolution = 2 * n^2 = {num_int_ops_direct} ops")
print(f"Time for convolution = {num_int_ops_direct} * {time_int_op} ns/op = {time_for_convolution} ns\n")
print(f"Total Time (Direct) = Time for conversion + Time for convolution")
print(f"Total Time (Direct) = {time_for_conversion} ns + {time_for_convolution} ns = {T_direct} ns\n")

# Step 4: Compare the two times and conclude
print("--- Comparison and Conclusion ---")
print(f"FFT-based Algorithm Time: {T_fft} ns")
print(f"Direct Integer Algorithm Time: {T_direct} ns")

if T_fft < T_direct:
    print("The original FFT-based algorithm is faster.")
    answer = "Y"
else:
    print("The direct integer algorithm is faster.")
    answer = "N"

print(f"\nIs the original FFT-based algorithm already the fastest? (Y/N)")
print(f"<<<{answer}>>>")