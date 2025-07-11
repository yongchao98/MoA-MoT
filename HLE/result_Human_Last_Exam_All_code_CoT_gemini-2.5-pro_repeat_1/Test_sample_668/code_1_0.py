import math

def solve_convolution_time():
    """
    Calculates and compares the time for different convolution algorithms.
    """
    # Machine operation times in nanoseconds (ns)
    int_add_time = 1
    int_mul_time = 2
    fp_add_time = 9
    fp_mul_time = 19

    # Size of the series
    n = 1000

    print("To determine the fastest algorithm, we will estimate the total execution time for three different approaches:")
    print("A. FFT-based convolution with floating points")
    print("B. Direct convolution with integers")
    print("C. Direct convolution with floating points")
    print("-" * 50)

    # --- Step 1: Direct Convolution ---
    print("\nStep 1: Calculate the time for Direct Convolution.")
    print(f"For two series of length n={n}, direct convolution requires n^2 multiplications and n*(n-1) additions.")
    
    num_mult_direct = n * n
    num_add_direct = n * (n - 1)

    print(f"Number of multiplications = {n}^2 = {num_mult_direct:,}")
    print(f"Number of additions = {n} * ({n} - 1) = {num_add_direct:,}")
    
    # Case B: Direct convolution with integers
    print("\nAnalysis for B: Direct convolution with integers")
    time_mult_direct_int = num_mult_direct * int_mul_time
    time_add_direct_int = num_add_direct * int_add_time
    total_time_direct_int = time_mult_direct_int + time_add_direct_int
    print(f"Time = ({num_mult_direct:,} multiplications * {int_mul_time} ns/mult) + ({num_add_direct:,} additions * {int_add_time} ns/add)")
    print(f"Time = {time_mult_direct_int:,} ns + {time_add_direct_int:,} ns = {total_time_direct_int:,} ns")

    # Case C: Direct convolution with floating points
    print("\nAnalysis for C: Direct convolution with floating points")
    time_mult_direct_fp = num_mult_direct * fp_mul_time
    time_add_direct_fp = num_add_direct * fp_add_time
    total_time_direct_fp = time_mult_direct_fp + time_add_direct_fp
    print(f"Time = ({num_mult_direct:,} multiplications * {fp_mul_time} ns/mult) + ({num_add_direct:,} additions * {fp_add_time} ns/add)")
    print(f"Time = {time_mult_direct_fp:,} ns + {time_add_direct_fp:,} ns = {total_time_direct_fp:,} ns")
    print("-" * 50)

    # --- Step 2: FFT-based Convolution ---
    print("\nStep 2: Calculate the time for FFT-based Convolution.")
    print("This method involves 2 forward FFTs, 1 element-wise complex multiplication, and 1 inverse FFT.")
    
    conv_len = 2 * n - 1
    N = 2**math.ceil(math.log2(conv_len))
    log2N = int(math.log2(N))

    print(f"The length of the convolution is 2*n - 1 = {conv_len}. We use an FFT size N that is the next power of 2.")
    print(f"So, N = {N} and log2(N) = {log2N}.")
    
    print("\nAn N-point complex FFT requires (N/2)*log2(N) complex multiplications and N*log2(N) complex additions.")
    print("A complex multiplication takes 4 real multiplications and 2 real additions.")
    print("A complex addition takes 2 real additions.")

    # Operations for one N-point FFT
    real_mult_per_fft = (N / 2) * log2N * 4
    real_add_per_fft = ((N / 2) * log2N * 2) + (N * log2N * 2)
    
    # Operations for element-wise complex product
    elementwise_mult = N * 4
    elementwise_add = N * 2

    # Total operations for the whole FFT-based convolution (3 FFTs + 1 element-wise product)
    total_fft_mult = 3 * real_mult_per_fft + elementwise_mult
    total_fft_add = 3 * real_add_per_fft + elementwise_add
    
    print("\nTotal operations for the entire process:")
    print(f"Total real multiplications = 3 * (({N}/2)*{log2N}*4) + ({N}*4) = {int(total_fft_mult):,}")
    print(f"Total real additions = 3 * (({N}/2)*{log2N}*2 + {N}*{log2N}*2) + ({N}*2) = {int(total_fft_add):,}")

    # Case A: FFT-based convolution with floating points
    print("\nAnalysis for A: FFT-based convolution with floating points")
    time_mult_fft_fp = total_fft_mult * fp_mul_time
    time_add_fft_fp = total_fft_add * fp_add_time
    total_time_fft_fp = time_mult_fft_fp + time_add_fft_fp
    print(f"Time = ({int(total_fft_mult):,} multiplications * {fp_mul_time} ns/mult) + ({int(total_fft_add):,} additions * {fp_add_time} ns/add)")
    print(f"Time = {int(time_mult_fft_fp):,} ns + {int(time_add_fft_fp):,} ns = {int(total_time_fft_fp):,} ns")
    print("-" * 50)

    # --- Step 3: Comparison ---
    print("\nStep 3: Compare the results.")
    print(f"Time for Direct Convolution (Integer):      {total_time_direct_int:15,} ns")
    print(f"Time for FFT-based Convolution (Float):     {int(total_time_fft_fp):15,} ns")
    print(f"Time for Direct Convolution (Float):      {total_time_direct_fp:15,} ns")

    print("\nConclusion: Direct convolution with integers is the fastest method for n=1000 on this machine.")

if __name__ == "__main__":
    solve_convolution_time()