import math

# --- Given parameters ---
n = 1000
T_float = 5  # ns, time for a floating point operation
T_int = 1    # ns, time for an integer operation
T_call = 15  # ns, time to call a function

# --- Calculation for FFT-based algorithm ---
# The total time is derived from the recurrence T(n) = 2*T(n/2) + (4*n*T_float + T_call)
# which solves to (4*n*log2(n))*T_float + (n-1)*T_call
log2_n = math.log2(n)
fft_time = (4 * n * log2_n) * T_float + (n - 1) * T_call

# --- Calculation for Integer-based algorithm ---
# The time is the sum of conversion cost and convolution cost
int_time = (2 * n) * T_float + (2 * n**2) * T_int

# --- Output the analysis and results ---
print("--- Algorithm Performance Analysis ---")
print(f"Parameters: n={n}, T_float={T_float} ns, T_int={T_int} ns, T_call={T_call} ns\n")

print("1. FFT-based Algorithm Cost")
print("   Formula: (4 * n * log2(n)) * T_float + (n - 1) * T_call")
# Print the equation with plugged-in numbers
print(f"   Time = (4 * {n} * {log2_n:.4f}) * {T_float} + ({n - 1}) * {T_call}")
# Print the value of each term before final summation
term1_val = 4 * n * log2_n * T_float
term2_val = (n - 1) * T_call
print(f"   Time = {term1_val:.2f} ns + {term2_val:.2f} ns")
print(f"   Total Time (FFT) = {fft_time:.2f} ns\n")


print("2. Integer-based Algorithm Cost")
print("   Formula: (2 * n) * T_float + (2 * n^2) * T_int")
# Print the equation with plugged-in numbers
print(f"   Time = (2 * {n}) * {T_float} + (2 * {n}**2) * {T_int}")
# Print the value of each term before final summation
term3_val = (2 * n) * T_float
term4_val = (2 * n**2) * T_int
print(f"   Time = {term3_val:.2f} ns + {term4_val:.2f} ns")
print(f"   Total Time (Integer) = {int_time:.2f} ns\n")

print("--- Conclusion ---")
print(f"FFT Algorithm Time: {fft_time:.2f} ns")
print(f"Integer Algorithm Time: {int_time:.2f} ns")

if fft_time < int_time:
    result = "Y"
    print("The original FFT-based algorithm is faster.")
else:
    result = "N"
    print("The proposed integer-based algorithm is faster.")

print(f"<<<{result}>>>")