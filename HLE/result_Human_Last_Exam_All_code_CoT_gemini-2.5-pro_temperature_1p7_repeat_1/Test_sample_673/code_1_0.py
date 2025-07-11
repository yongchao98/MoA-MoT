import math

# Given parameters
n = 1000
float_op_time = 5  # in nanoseconds
int_op_time = 1    # in nanoseconds
func_call_time = 15 # in nanoseconds

# --- FFT-based Algorithm ---
# A divide-and-conquer algorithm like FFT involves a recursion tree.
# For n inputs, this results in approximately 2n - 1 function calls.
fft_num_calls = 2 * n - 1
fft_num_float_ops = 4 * n
fft_cost_calls = fft_num_calls * func_call_time
fft_cost_ops = fft_num_float_ops * float_op_time
total_time_fft = fft_cost_calls + fft_cost_ops

# --- Direct Integer Convolution Algorithm ---
# This is a non-recursive approach, so we assume 1 function call.
direct_num_calls = 1
direct_num_float_ops = 2 * n
direct_num_int_ops = 2 * n**2
direct_cost_calls = direct_num_calls * func_call_time
direct_cost_float_ops = direct_num_float_ops * float_op_time
direct_cost_int_ops = direct_num_int_ops * int_op_time
total_time_direct = direct_cost_calls + direct_cost_float_ops + direct_cost_int_ops

print("### Algorithm Cost Analysis ###")
print(f"Vector size (n) = {n}\n")

# --- Print FFT Time Calculation ---
print("--- FFT-based Algorithm Time Calculation ---")
print("T_fft = (2 * n - 1) * func_call_time + (4 * n) * float_op_time")
print(f"T_fft = (2 * {n} - 1) * {func_call_time} ns + (4 * {n}) * {float_op_time} ns")
print(f"T_fft = ({fft_num_calls}) * {func_call_time} ns + ({fft_num_float_ops}) * {float_op_time} ns")
print(f"T_fft = {fft_cost_calls} ns + {fft_cost_ops} ns")
print(f"Total FFT Time = {total_time_fft} ns\n")

# --- Print Direct Convolution Time Calculation ---
print("--- Direct Integer Convolution Time Calculation ---")
print("T_direct = 1 * func_call_time + (2 * n) * float_op_time + (2 * n^2) * int_op_time")
print(f"T_direct = {direct_num_calls} * {func_call_time} ns + (2 * {n}) * {float_op_time} ns + (2 * {n}*{n}) * {int_op_time} ns")
print(f"T_direct = {direct_cost_calls} ns + ({direct_num_float_ops}) * {float_op_time} ns + ({direct_num_int_ops}) * {int_op_time} ns")
print(f"T_direct = {direct_cost_calls} ns + {direct_cost_float_ops} ns + {direct_cost_int_ops} ns")
print(f"Total Direct Time = {total_time_direct} ns\n")

# --- Conclusion ---
print("--- Conclusion ---")
if total_time_fft < total_time_direct:
    print(f"FFT-based time ({total_time_fft} ns) is less than direct convolution time ({total_time_direct} ns).")
    answer = "Y"
else:
    print(f"Direct convolution time ({total_time_direct} ns) is less than FFT-based time ({total_time_fft} ns).")
    answer = "N"

print("\nIs the original FFT-based algorithm is already the fastest?")
print(f"<<<{answer}>>>")