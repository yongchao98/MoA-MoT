import math

# Define the given parameters
n = 1000
t_fp = 5      # Time for a floating point operation in nanoseconds
t_int = 1     # Time for an integer operation in nanoseconds
t_call = 15   # Time for a function call in nanoseconds

# --- Calculation for the FFT-based Algorithm ---
print("--- FFT-based Algorithm Time Calculation ---")
print("The time complexity follows the recurrence T(n) = 2*T(n/2) + 4*n*t_fp + 2*t_call.")
print("The closed-form solution is T(n) = (4 * t_fp * n * log2(n)) + (2 * t_call * (n-1)).\n")

# Calculate the terms of the FFT time formula
fft_fp_ops_term = 4 * t_fp * n * math.log2(n)
fft_call_term = 2 * t_call * (n - 1)
total_time_fft = fft_fp_ops_term + fft_call_term

# Print the equation with numerical values
print(f"Time_FFT = (4 * {t_fp} * {n} * log2({n})) + (2 * {t_call} * ({n} - 1))")
print(f"Time_FFT = {fft_fp_ops_term:.2f} + {fft_call_term:.2f}")
print(f"Total Time for FFT-based algorithm = {total_time_fft:.2f} ns\n")


# --- Calculation for the Direct Integer-based Algorithm ---
print("--- Direct Integer-based Algorithm Time Calculation ---")
print("The time is the sum of costs for conversion, convolution, and one function call.")
print("Time = (2*n*t_fp) + (2*n^2*t_int) + t_call.\n")

# Calculate the terms of the direct algorithm time formula
direct_conversion_term = 2 * n * t_fp
direct_convolution_term = 2 * (n**2) * t_int
direct_call_term = t_call
total_time_direct = direct_conversion_term + direct_convolution_term + direct_call_term

# Print the equation with numerical values
print(f"Time_Direct = (2 * {n} * {t_fp}) + (2 * {n}**2 * {t_int}) + {t_call}")
print(f"Time_Direct = {direct_conversion_term} + {direct_convolution_term} + {direct_call_term}")
print(f"Total Time for Direct algorithm = {total_time_direct:.2f} ns\n")


# --- Comparison and Final Answer ---
print("--- Comparison ---")
print(f"FFT-based algorithm time: {total_time_fft:.2f} ns")
print(f"Direct algorithm time:      {total_time_direct:.2f} ns")

if total_time_fft < total_time_direct:
    print("The original FFT-based algorithm is faster.")
    answer = "Y"
else:
    print("The direct integer-based algorithm is faster or has equal performance.")
    answer = "N"

print("\nIs the original FFT-based algorithm is already the fastest?")
print(f"<<<{answer}>>>")