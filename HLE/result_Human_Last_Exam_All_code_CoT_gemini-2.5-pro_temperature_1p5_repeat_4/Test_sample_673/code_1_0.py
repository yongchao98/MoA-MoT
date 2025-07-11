import math

# Step 1: Define the given parameters
n = 1000  # vector size
t_flop = 5  # ns, time for a floating-point operation
t_int = 1   # ns, time for an integer operation
t_call = 15 # ns, time to call a function

print(f"Analysis for vector size n = {n}")
print("=" * 40)

# Step 2: Calculate the total time for the FFT-based algorithm (T1)
print("Algorithm 1: FFT-based Method")

# Cost from the "divide-and-conquer step", modeled as recursive function calls.
# A divide-and-conquer algorithm on size n (like recursive FFT) makes 2n-1 calls.
num_calls_fft = 2 * n - 1
cost_calls_fft = num_calls_fft * t_call
print(f"Cost of divide-and-conquer step (function calls):")
print(f"  {num_calls_fft} calls * {t_call} ns/call = {cost_calls_fft} ns")

# Cost from the floating point operations
num_flops_fft = 4 * n
cost_flops_fft = num_flops_fft * t_flop
print(f"Cost of floating-point operations:")
print(f"  {num_flops_fft} ops * {t_flop} ns/op = {cost_flops_fft} ns")

# Total time for Algorithm 1
total_time_fft = cost_calls_fft + cost_flops_fft
print(f"Total Time (T1) = {cost_calls_fft} + {cost_flops_fft} = {total_time_fft} ns")
print("-" * 40)


# Step 3: Calculate the total time for the direct integer algorithm (T2)
print("Algorithm 2: Direct Integer Convolution Method")

# This is an iterative algorithm, so we assume it is wrapped in a single function call.
num_calls_direct = 1
cost_calls_direct = num_calls_direct * t_call
print(f"Cost of function call:")
print(f"  {num_calls_direct} call * {t_call} ns/call = {cost_calls_direct} ns")

# Cost of converting real values to integers and back
num_flops_direct = 2 * n
cost_flops_direct = num_flops_direct * t_flop
print(f"Cost of floating-point conversions:")
print(f"  {num_flops_direct} ops * {t_flop} ns/op = {cost_flops_direct} ns")


# Cost of the direct convolution using integer operations
num_iops_direct = 2 * n**2
cost_iops_direct = num_iops_direct * t_int
print(f"Cost of integer convolution:")
print(f"  {num_iops_direct} ops * {t_int} ns/op = {cost_iops_direct} ns")


# Total time for Algorithm 2
total_time_direct = cost_calls_direct + cost_flops_direct + cost_iops_direct
print(f"Total Time (T2) = {cost_calls_direct} + {cost_flops_direct} + {cost_iops_direct} = {total_time_direct} ns")
print("=" * 40)


# Step 4 & 5: Compare the total times and conclude
print("Comparison:")
print(f"Time for FFT-based method (T1): {total_time_fft:,.0f} ns")
print(f"Time for direct integer method (T2): {total_time_direct:,.0f} ns")

if total_time_fft < total_time_direct:
    print("\nThe original FFT-based algorithm is faster.")
    answer = "Y"
else:
    print("\nThe direct integer convolution algorithm is faster.")
    answer = "N"

print("\nIs the original FFT-based algorithm is already the fastest?")
print(f"<<<{answer}>>>")