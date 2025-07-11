import math

# Step 1: Define the given parameters
n = 1000  # vector size
t_float = 5  # time for a floating-point operation in ns
t_int = 1    # time for an integer operation in ns
t_call = 15   # time for a function call in ns

print("Comparing two convolution algorithms based on the provided model.")
print(f"Parameters: n={n}, t_float={t_float} ns, t_int={t_int} ns, t_call={t_call} ns\n")

# Step 2: Calculate the execution time for the original FFT-based algorithm
# Model: Time = (2n - 1) * t_call + 4n * t_float
fft_calls = 2 * n - 1
fft_ops = 4 * n
time_fft = fft_calls * t_call + fft_ops * t_float

print("1. Calculating time for the original FFT-based algorithm:")
print(f"   Time = (Number of function calls * Time per call) + (Number of float ops * Time per float op)")
print(f"   Equation: (2 * {n} - 1) * {t_call} + (4 * {n}) * {t_float}")
print(f"   Calculation: ({fft_calls}) * {t_call} + ({fft_ops}) * {t_float} = {time_fft} ns\n")

# Step 3: Calculate the execution time for the direct integer-based algorithm
# Model: Time = 2n * t_float + 2n^2 * t_int
direct_conversion_ops = 2 * n
direct_convolution_ops = 2 * n**2
time_direct = direct_conversion_ops * t_float + direct_convolution_ops * t_int

print("2. Calculating time for the direct integer-based algorithm:")
print(f"   Time = (Number of conversion ops * Time per float op) + (Number of convolution ops * Time per int op)")
print(f"   Equation: (2 * {n}) * {t_float} + (2 * {n}^2) * {t_int}")
print(f"   Calculation: ({direct_conversion_ops}) * {t_float} + ({direct_convolution_ops}) * {t_int} = {time_direct} ns\n")

# Step 4: Compare the two results and provide the answer
print("3. Conclusion:")
print(f"   FFT-based algorithm time: {time_fft} ns")
print(f"   Direct integer algorithm time: {time_direct} ns")

if time_fft < time_direct:
    print(f"\nSince {time_fft} < {time_direct}, the original FFT-based algorithm is faster.")
    answer = "Y"
else:
    print(f"\nSince {time_fft} >= {time_direct}, the direct integer algorithm is faster or equal.")
    answer = "N"

print("\nIs the original FFT-based algorithm is already the fastest?")
print(f"<<<{answer}>>>")