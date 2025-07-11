import math

# --- Given Parameters ---
n = 1000  # Vector size
t_fp = 5  # Time for a floating-point operation in ns
t_int = 1 # Time for an integer operation in ns
t_call = 15 # Time for a function call in ns

# --- Algorithm 1: FFT-based ---
# Number of floating point operations
n_fp_fft = 4 * n
# Number of function calls for a divide-and-conquer algorithm splitting in two
n_calls_fft = 2 * n - 1
# Total time for the FFT-based algorithm
time_fft = n_fp_fft * t_fp + n_calls_fft * t_call

print("--- 1. FFT-based Algorithm Calculation ---")
print("The total time is the sum of time for floating-point operations and function calls.")
print(f"Time_FFT = (4 * n) * t_fp + (2 * n - 1) * t_call")
print(f"Time_FFT = (4 * {n}) * {t_fp} + (2 * {n} - 1) * {t_call}")
print(f"Time_FFT = {n_fp_fft} * {t_fp} + {n_calls_fft} * {t_call}")
print(f"Time_FFT = {n_fp_fft * t_fp} + {n_calls_fft * t_call}")
print(f"Time_FFT = {time_fft} ns\n")


# --- Algorithm 2: Integer-based Direct Convolution ---
# Number of floating point operations for conversions
n_fp_direct = 2 * n
# Number of integer operations for direct convolution
n_int_direct = 2 * n**2
# Total time for the direct integer algorithm
time_direct = n_fp_direct * t_fp + n_int_direct * t_int

print("--- 2. Integer-based Direct Algorithm Calculation ---")
print("The total time is the sum of time for floating-point conversions and integer operations.")
print(f"Time_Direct = (2 * n) * t_fp + (2 * n^2) * t_int")
print(f"Time_Direct = (2 * {n}) * {t_fp} + (2 * {n}**2) * {t_int}")
print(f"Time_Direct = {n_fp_direct} * {t_fp} + {n_int_direct} * {t_int}")
print(f"Time_Direct = {n_fp_direct * t_fp} + {n_int_direct * t_int}")
print(f"Time_Direct = {time_direct} ns\n")


# --- 3. Comparison and Conclusion ---
print("--- 3. Conclusion ---")
print(f"FFT-based algorithm total time: {time_fft} ns")
print(f"Integer-based direct algorithm total time: {time_direct} ns")

if time_fft < time_direct:
    is_fastest = "Y"
    print("\nThe original FFT-based algorithm is faster.")
else:
    is_fastest = "N"
    print("\nThe integer-based direct algorithm is faster or has the same speed.")

print(f"<<<{is_fastest}>>>")