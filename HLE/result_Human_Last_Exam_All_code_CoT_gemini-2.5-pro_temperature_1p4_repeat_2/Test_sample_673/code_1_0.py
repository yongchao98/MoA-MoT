# Define the given parameters
n = 1000  # vector size
t_float = 5  # ns per floating point operation
t_int = 1    # ns per integer operation
t_call = 15  # ns per function call

# --- Calculate the time for the FFT-based algorithm ---
print("Calculating time for FFT-based algorithm:")

# Number of function calls for a divide-and-conquer FFT is proportional to n.
# A common recurrence C(n) = 2C(n/2) + 1 leads to 2n-1 calls. We use 2n as an estimate.
calls_fft = 2 * n
time_calls_fft = calls_fft * t_call
print(f"Number of function calls: 2 * {n} = {calls_fft}")
print(f"Time for function calls: {calls_fft} * {t_call} ns = {time_calls_fft} ns")

# Number of floating point operations
ops_fft = 4 * n
time_ops_fft = ops_fft * t_float
print(f"Number of floating point operations: 4 * {n} = {ops_fft}")
print(f"Time for floating point operations: {ops_fft} * {t_float} ns = {time_ops_fft} ns")

# Total time for FFT-based algorithm
total_time_fft = time_calls_fft + time_ops_fft
print(f"Total time for FFT-based algorithm: {time_calls_fft} + {time_ops_fft} = {total_time_fft} ns")
print("-" * 30)

# --- Calculate the time for the direct integer-based algorithm ---
print("Calculating time for direct integer-based algorithm:")

# Number of function calls (can be implemented in one function)
calls_direct = 1
time_calls_direct = calls_direct * t_call
print(f"Number of function calls: {calls_direct}")
print(f"Time for function calls: {calls_direct} * {t_call} ns = {time_calls_direct} ns")

# Number of floating point operations for conversion
float_ops_direct = 2 * n
time_float_ops_direct = float_ops_direct * t_float
print(f"Number of floating point operations for conversion: 2 * {n} = {float_ops_direct}")
print(f"Time for floating point operations: {float_ops_direct} * {t_float} ns = {time_float_ops_direct} ns")

# Number of integer operations for convolution
int_ops_direct = 2 * n**2
time_int_ops_direct = int_ops_direct * t_int
print(f"Number of integer operations for convolution: 2 * {n}^2 = {int_ops_direct}")
print(f"Time for integer operations: {int_ops_direct} * {t_int} ns = {time_int_ops_direct} ns")

# Total time for direct integer-based algorithm
total_time_direct = time_calls_direct + time_float_ops_direct + time_int_ops_direct
print(f"Total time for direct integer-based algorithm: {time_calls_direct} + {time_float_ops_direct} + {time_int_ops_direct} = {total_time_direct} ns")
print("-" * 30)

# --- Compare the two algorithms ---
print("Comparison:")
print(f"FFT-based algorithm time: {total_time_fft} ns")
print(f"Direct integer-based algorithm time: {total_time_direct} ns")

if total_time_fft < total_time_direct:
    print("\nThe original FFT-based algorithm is faster.")
    answer = "Y"
else:
    print("\nThe proposed direct integer-based algorithm is faster.")
    answer = "N"

print(f"\nIs the original FFT-based algorithm is already the fastest? {answer}")
print(f"\n<<<Y>>>")