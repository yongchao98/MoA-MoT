import math

# Step 1: Define Constants from the problem description
n = 1000      # typical vector size
T_fp = 5      # ns, time for a floating point operation
T_int = 1     # ns, time for an integer operation
T_call = 15   # ns, time to call a function

# Step 2: Calculate the total time for the FFT-based algorithm
print("--- Calculating Time for FFT-based Algorithm ---")
print("Formula: T_fft = (2 * n - 1) * T_call + 4 * n * T_fp")

# Breaking down the calculation
fft_calls_term_val = (2 * n - 1) * T_call
fft_fp_term_val = 4 * n * T_fp
t_fft = fft_calls_term_val + fft_fp_term_val

# Printing the equation with substituted values
print(f"T_fft = (2 * {n} - 1) * {T_call} + 4 * {n} * {T_fp}")
print(f"T_fft = ({2*n-1}) * {T_call} + {4*n} * {T_fp}")
print(f"T_fft = {fft_calls_term_val} + {fft_fp_term_val}")
print(f"T_fft = {t_fft} ns\n")


# Step 3: Calculate the total time for the direct integer-based algorithm
print("--- Calculating Time for Direct Integer-based Algorithm ---")
print("Formula: T_direct = 2 * n * T_fp + 2 * n^2 * T_int")

# Breaking down the calculation
direct_fp_term_val = 2 * n * T_fp
n_squared = n**2
direct_int_term_val = 2 * n_squared * T_int
t_direct = direct_fp_term_val + direct_int_term_val

# Printing the equation with substituted values
print(f"T_direct = 2 * {n} * {T_fp} + 2 * {n}^2 * {T_int}")
print(f"T_direct = 2 * {n} * {T_fp} + 2 * {n_squared} * {T_int}")
print(f"T_direct = {direct_fp_term_val} + {direct_int_term_val}")
print(f"T_direct = {t_direct} ns\n")


# Step 4: Compare the results
print("--- Comparison ---")
if t_fft < t_direct:
    print(f"Time_FFT ({t_fft} ns) < Time_Direct ({t_direct} ns)")
    print("Conclusion: The original FFT-based algorithm is faster.")
    answer = 'Y'
else:
    print(f"Time_FFT ({t_fft} ns) >= Time_Direct ({t_direct} ns)")
    print("Conclusion: The direct integer-based algorithm is faster or has the same speed.")
    answer = 'N'

# Step 5: Final Answer
print("\nIs the original FFT-based algorithm is already the fastest?")
print(f"<<<{answer}>>>")