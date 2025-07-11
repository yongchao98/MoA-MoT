import math

# Step 1: Define Constants from the problem description
n = 1000
t_float_op = 5  # Time for one floating point operation in nanoseconds
t_int_op = 1    # Time for one integer operation in nanoseconds
t_call = 15     # Time for one function call in nanoseconds

print("--- Step 1: Provided Parameters ---")
print(f"Vector size (n): {n}")
print(f"Time per floating-point op: {t_float_op} ns")
print(f"Time per integer op: {t_int_op} ns")
print(f"Time per function call: {t_call} ns")
print("-" * 30 + "\n")


# Step 2: Calculate total time for Algorithm 1 (FFT-based)
print("--- Step 2: Algorithm 1 (FFT-based) Time Calculation ---")

# The divide-and-conquer step is modeled as function call overhead.
# A recursive FFT algorithm on a vector of size n has 2n-1 calls. We use 2n for approximation.
num_calls_algo1 = 2 * n
cost_calls_algo1 = num_calls_algo1 * t_call

# The problem states it performs 4n floating point operations.
num_float_ops_algo1 = 4 * n
cost_ops_algo1 = num_float_ops_algo1 * t_float_op

# Total time is the sum of the call overhead and operation cost.
total_time_algo1 = cost_calls_algo1 + cost_ops_algo1

print("Equation for total time: (2 * n * t_call) + (4 * n * t_float_op)")
print(f"Time for function calls = (2 * {n} * {t_call}) = {cost_calls_algo1} ns")
print(f"Time for float operations = (4 * {n} * {t_float_op}) = {cost_ops_algo1} ns")
print(f"Total time for Algorithm 1 = {cost_calls_algo1} + {cost_ops_algo1} = {total_time_algo1} ns")
print("-" * 30 + "\n")


# Step 3: Calculate total time for Algorithm 2 (Direct Integer Convolution)
print("--- Step 3: Algorithm 2 (Direct Integer Convolution) Time Calculation ---")

# Cost for converting real vectors to fixed-point integers
num_float_ops_algo2 = 2 * n
cost_conversion_algo2 = num_float_ops_algo2 * t_float_op

# Cost for direct convolution using integer operations
num_int_ops_algo2 = 2 * n**2
cost_convolution_algo2 = num_int_ops_algo2 * t_int_op

# Assume a single function call to perform the whole operation.
cost_call_algo2 = 1 * t_call

# Total time is the sum of conversion, convolution, and call overhead.
total_time_algo2 = cost_conversion_algo2 + cost_convolution_algo2 + cost_call_algo2

print("Equation for total time: (2 * n * t_float_op) + (2 * n^2 * t_int_op) + t_call")
print(f"Time for conversion = (2 * {n} * {t_float_op}) = {cost_conversion_algo2} ns")
print(f"Time for convolution = (2 * {n}^2 * {t_int_op}) = {cost_convolution_algo2} ns")
print(f"Time for single function call = {cost_call_algo2} ns")
print(f"Total time for Algorithm 2 = {cost_conversion_algo2} + {cost_convolution_algo2} + {cost_call_algo2} = {total_time_algo2} ns")
print("-" * 30 + "\n")


# Step 4: Compare the two algorithms
print("--- Step 4: Comparison ---")
print(f"Algorithm 1 (FFT-based) total time: {total_time_algo1} ns")
print(f"Algorithm 2 (Direct Integer) total time: {total_time_algo2} ns")

if total_time_algo1 < total_time_algo2:
    print("\nConclusion: The original FFT-based algorithm is faster.")
    answer = "Y"
else:
    print("\nConclusion: The direct integer convolution algorithm is faster.")
    answer = "N"

print("\nIs the original FFT-based algorithm is already the fastest?")
print(f"<<<{answer}>>>")