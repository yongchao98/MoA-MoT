import math

def calculate_convolution_times():
    """
    Calculates and compares the estimated execution time for convolution
    using different algorithms on a specified machine.
    """
    # Machine operation times in nanoseconds
    int_add_time = 1
    int_mul_time = 2
    fp_add_time = 9
    fp_mul_time = 19

    # Size of the input series
    n = 1000

    print("--- Analysis of Convolution Algorithms (n=1000) ---")
    
    # --- B. Direct Convolution with Integers ---
    print("\n--- Method B: Direct Convolution with Integers ---")
    direct_mults = n ** 2
    direct_adds = n * (n - 1)
    
    time_direct_int_mult = direct_mults * int_mul_time
    time_direct_int_add = direct_adds * int_add_time
    total_time_direct_int = time_direct_int_mult + time_direct_int_add
    
    print(f"Number of integer multiplications: {direct_mults}")
    print(f"Number of integer additions: {direct_adds}")
    print("Equation: Total time = (Number of multiplications * Integer multiplication time) + (Number of additions * Integer addition time)")
    print(f"Total time = {direct_mults} * {int_mul_time} ns + {direct_adds} * {int_add_time} ns = {total_time_direct_int:,} ns")

    # --- C. Direct Convolution with Floating Points ---
    print("\n--- Method C: Direct Convolution with Floating Points ---")
    time_direct_fp_mult = direct_mults * fp_mul_time
    time_direct_fp_add = direct_adds * fp_add_time
    total_time_direct_fp = time_direct_fp_mult + time_direct_fp_add
    
    print(f"Number of floating point multiplications: {direct_mults}")
    print(f"Number of floating point additions: {direct_adds}")
    print("Equation: Total time = (Number of multiplications * FP multiplication time) + (Number of additions * FP addition time)")
    print(f"Total time = {direct_mults} * {fp_mul_time} ns + {direct_adds} * {fp_add_time} ns = {total_time_direct_fp:,} ns")

    # --- A. FFT-based Convolution ---
    print("\n--- Method A: FFT-based Convolution (Floating Point) ---")
    fft_size_n = 2 * n - 1
    # Find the next power of 2
    N = 1 << (fft_size_n - 1).bit_length()
    log2_N = math.log2(N)

    # Operations for 3 FFTs/IFFTs and 1 pointwise multiplication
    # A complex multiplication is 4 real muls and 2 real adds
    # An FFT of size N has (N/2)*log2(N) complex muls and N*log2(N) complex adds
    
    # Total real multiplications: 3 * (N/2 * log2(N) * 4) + (N * 4)
    fft_total_mults = 3 * (N / 2 * log2_N * 4) + (N * 4)
    
    # Total real additions: 3 * ( (N/2*log2(N)*2) + (N*log2(N)*2) ) + (N * 2)
    fft_total_adds = 3 * ((N / 2 * log2_N * 2) + (N * log2_N * 2)) + (N * 2)
    
    # Cast to int for clean printing
    fft_total_mults = int(fft_total_mults)
    fft_total_adds = int(fft_total_adds)

    time_fft_mult = fft_total_mults * fp_mul_time
    time_fft_add = fft_total_adds * fp_add_time
    total_time_fft = time_fft_mult + time_fft_add

    print(f"FFT Size (N): {N}")
    print(f"Number of floating point multiplications: {fft_total_mults}")
    print(f"Number of floating point additions: {fft_total_adds}")
    print("Equation: Total time = (Number of multiplications * FP multiplication time) + (Number of additions * FP addition time)")
    print(f"Total time = {fft_total_mults} * {fp_mul_time} ns + {fft_total_adds} * {fp_add_time} ns = {total_time_fft:,} ns")
    
    print("\n--- Conclusion ---")
    print(f"Time (Direct Integer): {total_time_direct_int:,} ns")
    print(f"Time (FFT Float):      {total_time_fft:,} ns")
    print(f"Time (Direct Float):   {total_time_direct_fp:,} ns")
    print("\nThe fastest algorithm is Direct Convolution with Integers.")


if __name__ == '__main__':
    calculate_convolution_times()