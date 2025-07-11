import math

# --- Machine and Problem Parameters ---
n = 1000
int_add_time = 1  # ns
int_mult_time = 2  # ns
float_add_time = 9  # ns
float_mult_time = 19 # ns

print("Analysis of Convolution Algorithms for n = 1000")
print("=================================================")
print(f"Machine Operation Times (ns):")
print(f"- Integer Add/Sub: {int_add_time}")
print(f"- Integer Mult: {int_mult_time}")
print(f"- Floating Point Add/Sub: {float_add_time}")
print(f"- Floating Point Mult: {float_mult_time}")
print("\n---")


# --- Method 1: Direct Convolution ---
print("1. Direct Convolution Calculation")
num_mult_direct = n**2
num_add_direct = (n - 1)**2

print(f"   Number of multiplications: {n}^2 = {num_mult_direct}")
print(f"   Number of additions: ({n}-1)^2 = {num_add_direct}")

# B. Direct convolution with integers
time_direct_int = (num_mult_direct * int_mult_time) + (num_add_direct * int_add_time)
print(f"\n   Time with Integers:")
print(f"   ({num_mult_direct} * {int_mult_time} ns) + ({num_add_direct} * {int_add_time} ns) = {time_direct_int:,} ns")

# C. Direct convolution with floating points
time_direct_float = (num_mult_direct * float_mult_time) + (num_add_direct * float_add_time)
print(f"\n   Time with Floating Points:")
print(f"   ({num_mult_direct} * {float_mult_time} ns) + ({num_add_direct} * {float_add_time} ns) = {time_direct_float:,} ns")
print("\n---")


# --- Method 2: FFT-based Convolution ---
print("2. FFT-based Convolution Calculation")
# Find the next power of 2 for FFT size N >= 2n-1
fft_len_min = 2 * n - 1
N = 2**math.ceil(math.log2(fft_len_min))
log2_N = math.log2(N)

print(f"   Input series length n = {n}")
print(f"   Required padded length >= 2*n-1 = {fft_len_min}")
print(f"   FFT size N (next power of 2) = {N}")
print(f"   log2(N) = {int(log2_N)}")

# Operations for FFT/IFFT
# A single complex multiplication (a+ib)(c+id) = (ac-bd) + i(ad+bc) takes:
# 4 real multiplications and 2 real additions.
# A single complex addition (a+ib)+(c+id) = (a+c)+i(b+d) takes:
# 2 real additions.
real_mult_per_complex_mult = 4
real_add_per_complex_mult = 2
real_add_per_complex_add = 2

# Total operations = 2 FFTs + 1 IFFT + element-wise product
num_transforms = 3
complex_mult_in_transforms = num_transforms * (N / 2) * log2_N
complex_add_in_transforms = num_transforms * N * log2_N
complex_mult_elementwise = N

total_complex_mult = complex_mult_in_transforms + complex_mult_elementwise
total_complex_add = complex_add_in_transforms

# Convert to total real floating point operations
total_float_mult = total_complex_mult * real_mult_per_complex_mult
total_float_add = (total_complex_mult * real_add_per_complex_mult) + \
                  (total_complex_add * real_add_per_complex_add)

print(f"\n   Total complex multiplications (3 transforms + element-wise): {int(total_complex_mult)}")
print(f"   Total complex additions (3 transforms): {int(total_complex_add)}")

print(f"\n   Equivalent real floating-point operations:")
print(f"   Total multiplications = {int(total_complex_mult)} * {real_mult_per_complex_mult} = {int(total_float_mult)}")
print(f"   Total additions = ({int(total_complex_mult)} * {real_add_per_complex_mult}) + ({int(total_complex_add)} * {real_add_per_complex_add}) = {int(total_float_add)}")

# A. FFT time calculation
time_fft = (total_float_mult * float_mult_time) + (total_float_add * float_add_time)
print(f"\n   Time with FFT (Floating Point):")
print(f"   ({int(total_float_mult)} * {float_mult_time} ns) + ({int(total_float_add)} * {float_add_time} ns) = {time_fft:,.0f} ns")
print("\n---")


# --- Conclusion ---
print("3. Conclusion")
times = {
    "Direct convolution with integers": time_direct_int,
    "FFT-based convolution": time_fft,
    "Direct convolution with floating points": time_direct_float
}

fastest_method = min(times, key=times.get)
fastest_time = times[fastest_method]

print(f"   Time Comparison:")
print(f"   - Direct Integer: {time_direct_int:,} ns")
print(f"   - FFT (Float):    {time_fft:,.0f} ns")
print(f"   - Direct Float:   {time_direct_float:,} ns")
print(f"\n   The fastest method is '{fastest_method}' with an estimated time of {fastest_time:,.0f} ns.")
