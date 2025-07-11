import math

def solve_convolution_time():
    # Machine operation times in nanoseconds (ns)
    time_int_add = 1
    time_int_mult = 2
    time_float_add = 9
    time_float_mult = 19

    # Series length
    n = 1000

    print("Step-by-step calculation to find the fastest convolution algorithm for n=1000.")
    print("-" * 70)

    # --- Case 1: Direct Convolution with Integers ---
    print("Analysis for Direct Convolution with Integers:")
    n_ops_mult_direct = n**2
    n_ops_add_direct = (n - 1)**2
    print(f"Direct convolution requires n^2 multiplications and (n-1)^2 additions.")
    print(f"For n = {n}:")
    print(f"  Number of integer multiplications = {n}^2 = {n_ops_mult_direct}")
    print(f"  Number of integer additions = ({n}-1)^2 = {n_ops_add_direct}")

    time_direct_int_mult = n_ops_mult_direct * time_int_mult
    time_direct_int_add = n_ops_add_direct * time_int_add
    total_time_direct_int = time_direct_int_mult + time_direct_int_add

    print(f"Total time = (Number of multiplications * time) + (Number of additions * time)")
    print(f"Total time = {n_ops_mult_direct} * {time_int_mult} ns + {n_ops_add_direct} * {time_int_add} ns")
    print(f"Total time = {time_direct_int_mult} ns + {time_direct_int_add} ns = {total_time_direct_int} ns")
    print("-" * 70)

    # --- Case 2: Direct Convolution with Floating Points ---
    print("Analysis for Direct Convolution with Floating Points:")
    # The number of operations is the same as the integer case.
    print(f"The number of operations is the same as for integers:")
    print(f"  Number of floating point multiplications = {n_ops_mult_direct}")
    print(f"  Number of floating point additions = {n_ops_add_direct}")
    
    time_direct_float_mult = n_ops_mult_direct * time_float_mult
    time_direct_float_add = n_ops_add_direct * time_float_add
    total_time_direct_float = time_direct_float_mult + time_direct_float_add
    
    print(f"Total time = (Number of multiplications * time) + (Number of additions * time)")
    print(f"Total time = {n_ops_mult_direct} * {time_float_mult} ns + {n_ops_add_direct} * {time_float_add} ns")
    print(f"Total time = {time_direct_float_mult} ns + {time_direct_float_add} ns = {total_time_direct_float} ns")
    print("-" * 70)

    # --- Case 3: FFT-based Convolution ---
    print("Analysis for FFT-based Convolution (with Floating Points):")
    conv_len = 2 * n - 1
    # Find the next power of 2 >= conv_len
    N = 1
    while N < conv_len:
        N *= 2
    log2_N = int(math.log2(N))

    print(f"The length of the convolved sequence is 2*n - 1 = {conv_len}.")
    print(f"Data is padded to the next power of 2, so the FFT size N = {N}.")
    print(f"The calculation involves 2 forward FFTs, 1 pointwise multiplication, and 1 inverse FFT (IFFT).")

    # Operations for one N-point FFT
    # A complex multiplication (a+ib)(c+id) requires 4 real multiplications and 2 real additions.
    # A complex addition (a+ib)+(c+id) requires 2 real additions.
    complex_mults_per_fft = (N / 2) * log2_N
    complex_adds_per_fft = N * log2_N
    
    real_mults_per_fft = int(complex_mults_per_fft * 4)
    real_adds_per_fft = int(complex_mults_per_fft * 2 + complex_adds_per_fft * 2)

    print(f"\nOperations for one {N}-point FFT:")
    print(f"  Number of real multiplications = ({N}/2 * log2({N})) * 4 = {real_mults_per_fft}")
    print(f"  Number of real additions = (({N}/2 * log2({N}))*2) + (N*log2({N})*2) = {real_adds_per_fft}")

    # Total operations for the whole process
    # 2 FFTs, 1 IFFT = 3 FFT-like operations
    # N pointwise complex multiplications
    total_fft_ops_mults = 3 * real_mults_per_fft + N * 4
    total_fft_ops_adds = 3 * real_adds_per_fft + N * 2

    print(f"\nTotal operations for the complete algorithm (2 FFTs, 1 pointwise mult, 1 IFFT):")
    print(f"  Total float multiplications = 3 * {real_mults_per_fft} + {N} * 4 = {total_fft_ops_mults}")
    print(f"  Total float additions = 3 * {real_adds_per_fft} + {N} * 2 = {total_fft_ops_adds}")

    time_fft_mult = total_fft_ops_mults * time_float_mult
    time_fft_add = total_fft_ops_adds * time_float_add
    total_time_fft = time_fft_mult + time_fft_add
    
    print(f"\nTotal time = (Total multiplications * time) + (Total additions * time)")
    print(f"Total time = {total_fft_ops_mults} * {time_float_mult} ns + {total_fft_ops_adds} * {time_float_add} ns")
    print(f"Total time = {time_fft_mult} ns + {time_fft_add} ns = {total_time_fft} ns")
    print("-" * 70)

    # --- Comparison and Conclusion ---
    print("Comparison of Total Times:")
    print(f"  Direct convolution with integers:       {total_time_direct_int} ns")
    print(f"  Direct convolution with floating points:  {total_time_direct_float} ns")
    print(f"  FFT-based convolution:                  {total_time_fft} ns")
    print("-" * 70)
    
    times = {
        "Direct convolution with integers": total_time_direct_int,
        "Direct convolution with floating points": total_time_direct_float,
        "FFT": total_time_fft,
    }

    fastest_method_name = min(times, key=times.get)
    print(f"Conclusion: The fastest algorithm is '{fastest_method_name}'.")

# Execute the function to print the analysis
solve_convolution_time()