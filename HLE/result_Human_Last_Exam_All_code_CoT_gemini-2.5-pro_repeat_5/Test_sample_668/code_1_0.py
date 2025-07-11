import math

def calculate_convolution_times():
    """
    Calculates and compares the execution times for different convolution algorithms.
    """
    # Machine operation times in nanoseconds (ns)
    int_add_time = 1
    int_mul_time = 2
    fp_add_time = 9
    fp_mul_time = 19

    # Series length
    n = 1000

    print("--- Analysis for n = 1000 ---")

    # --- Direct Convolution ---
    print("\n1. Direct Convolution Calculation:")
    direct_mults = n * n
    direct_adds = (n - 1) * (n - 1)
    
    # Integer case
    time_direct_int = (direct_mults * int_mul_time) + (direct_adds * int_add_time)
    print("   - Integer-based:")
    print(f"     Equation: ({direct_mults} multiplications * {int_mul_time} ns) + ({direct_adds} additions * {int_add_time} ns)")
    print(f"     Time = {direct_mults * int_mul_time} ns + {direct_adds * int_add_time} ns = {int(time_direct_int)} ns")

    # Floating point case
    time_direct_fp = (direct_mults * fp_mul_time) + (direct_adds * fp_add_time)
    print("   - Floating point-based:")
    print(f"     Equation: ({direct_mults} multiplications * {fp_mul_time} ns) + ({direct_adds} additions * {fp_add_time} ns)")
    print(f"     Time = {direct_mults * fp_mul_time} ns + {direct_adds * fp_add_time} ns = {int(time_direct_fp)} ns")

    # --- FFT-based Convolution ---
    print("\n2. FFT-based Convolution Calculation:")
    
    # Find the next power of 2 for FFT size N
    N = 1
    while N < (2 * n - 1):
        N *= 2
    log2N = int(math.log2(N))

    print(f"   - FFT size N = {N} (next power of 2 >= 2*n-1 = {2*n-1})")

    # Operations for 3 FFTs (2 forward, 1 inverse)
    # Each complex mult = 4 real mults, 2 real adds
    # Each complex add = 2 real adds
    # Total real mults for 3 FFTs = 3 * (N/2 * log2N * 4) = 6 * N * log2N
    # Total real adds for 3 FFTs = 3 * (N/2 * log2N * 2 + N*log2N*2) = 9 * N * log2N
    fft_mults = 6 * N * log2N
    fft_adds = 9 * N * log2N
    
    # Operations for element-wise complex multiplication (N multiplications)
    pointwise_mults = N * 4
    pointwise_adds = N * 2

    total_fp_mults = fft_mults + pointwise_mults
    total_fp_adds = fft_adds + pointwise_adds
    
    time_fft = (total_fp_mults * fp_mul_time) + (total_fp_adds * fp_add_time)
    print("   - Floating point-based (FFT):")
    print(f"     Equation: ({total_fp_mults} multiplications * {fp_mul_time} ns) + ({total_fp_adds} additions * {fp_add_time} ns)")
    print(f"     Time = {total_fp_mults * fp_mul_time} ns + {total_fp_adds * fp_add_time} ns = {int(time_fft)} ns")

    # --- Conclusion ---
    print("\n--- Conclusion ---")
    times = {
        "Direct convolution with integers": time_direct_int,
        "Direct convolution with floating points": time_direct_fp,
        "FFT-based convolution": time_fft,
    }
    fastest_method = min(times, key=times.get)
    print(f"Comparing the total times:\n- Direct (Integer): {int(time_direct_int)} ns\n- Direct (Float):   {int(time_direct_fp)} ns\n- FFT (Float):      {int(time_fft)} ns")
    print(f"\nThe fastest algorithm is: {fastest_method}")


if __name__ == '__main__':
    calculate_convolution_times()