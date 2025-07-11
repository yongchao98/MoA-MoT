import math

# Step 1: Define the given parameters.
n = 1000
time_flop = 5  # ns per floating point operation
time_int_op = 1  # ns per integer operation
time_call = 15  # ns per function call

# Step 2: Calculate the total time for the FFT-based algorithm.
# The number of recursive function calls for a D&C algorithm like FFT is ~2n.
fft_calls = 2 * n
fft_flops = 4 * n

# Calculate the components of the total time for the FFT algorithm.
fft_call_time = fft_calls * time_call
fft_flop_time = fft_flops * time_flop
total_time_fft = fft_call_time + fft_flop_time

# Step 3: Calculate the total time for the integer-based direct convolution algorithm.
int_conv_flops = 2 * n
int_conv_ops = 2 * n**2

# Calculate the components of the total time for the integer algorithm.
int_conv_flop_time = int_conv_flops * time_flop
int_conv_op_time = int_conv_ops * time_int_op
total_time_int = int_conv_flop_time + int_conv_op_time

# Step 4: Print the calculations and the comparison.
print("--- Algorithm Performance Comparison ---")
print(f"\nParameters:\n- Vector size (n): {n}\n- FLOP time: {time_flop} ns\n- Integer op time: {time_int_op} ns\n- Function call time: {time_call} ns")

print("\n1. FFT-based Algorithm Time Calculation:")
print(f"Total Time = (Function Calls * Time per Call) + (FLOPS * Time per FLOP)")
print(f"Total Time = ({fft_calls} calls * {time_call} ns) + ({fft_flops} FLOPS * {time_flop} ns)")
print(f"Total Time = {fft_call_time} ns + {fft_flop_time} ns")
print(f"Total Time = {total_time_fft} ns")

print("\n2. Integer-based Direct Convolution Time Calculation:")
print(f"Total Time = (Conversion FLOPS * Time per FLOP) + (Integer Ops * Time per Integer Op)")
print(f"Total Time = ({int_conv_flops} FLOPS * {time_flop} ns) + ({int_conv_ops} ops * {time_int_op} ns)")
print(f"Total Time = {int_conv_flop_time} ns + {int_conv_op_time} ns")
print(f"Total Time = {total_time_int} ns")

print("\n--- Conclusion ---")
print(f"FFT-based Algorithm Time: {total_time_fft} ns")
print(f"Integer-based Algorithm Time: {total_time_int} ns")

is_fft_fastest = total_time_fft < total_time_int

if is_fft_fastest:
    print("\nThe original FFT-based algorithm is faster.")
    answer = "Y"
else:
    print("\nThe proposed integer-based algorithm is faster.")
    answer = "N"

print("\nIs the original FFT-based algorithm already the fastest?")
print(f"<<<{answer}>>>")