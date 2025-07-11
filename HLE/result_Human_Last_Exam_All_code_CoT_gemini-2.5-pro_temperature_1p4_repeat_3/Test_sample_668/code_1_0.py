import math

def calculate_convolution_times():
    # Given parameters
    n = 1000
    int_add_time = 1  # ns
    int_mult_time = 2  # ns
    float_add_time = 9  # ns
    float_mult_time = 19  # ns

    # --- Direct Convolution (Integers) ---
    print("--- Calculating Time for Direct Convolution (Integers) ---")
    num_mults_direct = n * n
    num_adds_direct = n * (n - 1)
    
    time_direct_int_mult = num_mults_direct * int_mult_time
    time_direct_int_add = num_adds_direct * int_add_time
    total_time_direct_int = time_direct_int_mult + time_direct_int_add
    
    print(f"Number of integer multiplications: {num_mults_direct}")
    print(f"Number of integer additions: {num_adds_direct}")
    print(f"Calculation: ({num_mults_direct} * {int_mult_time} ns) + ({num_adds_direct} * {int_add_time} ns)")
    print(f"Total time for Direct Convolution (Integers): {total_time_direct_int} ns\n")

    # --- Direct Convolution (Floating Points) ---
    print("--- Calculating Time for Direct Convolution (Floating Points) ---")
    time_direct_float_mult = num_mults_direct * float_mult_time
    time_direct_float_add = num_adds_direct * float_add_time
    total_time_direct_float = time_direct_float_mult + time_direct_float_add
    
    print(f"Number of floating point multiplications: {num_mults_direct}")
    print(f"Number of floating point additions: {num_adds_direct}")
    print(f"Calculation: ({num_mults_direct} * {float_mult_time} ns) + ({num_adds_direct} * {float_add_time} ns)")
    print(f"Total time for Direct Convolution (Floating Points): {total_time_direct_float} ns\n")

    # --- FFT-based Convolution (Floating Points) ---
    print("--- Calculating Time for FFT-based Convolution ---")
    # 1. Find N
    N = 1
    while N < (2 * n - 1):
        N *= 2
    log2_N = math.log2(N)
    
    # 2. Operations for one FFT
    # complex mult = 4 real mults + 2 real adds
    # complex add = 2 real adds
    num_complex_mults_fft = (N / 2) * log2_N
    num_complex_adds_fft = N * log2_N
    
    real_mults_per_fft = int(num_complex_mults_fft * 4)
    real_adds_per_fft = int((num_complex_mults_fft * 2) + (num_complex_adds_fft * 2))
    
    # 3. Total operations for convolution
    # 2 forward FFTs, 1 inverse FFT, 1 element-wise complex multiplication of size N
    total_real_mults_fft = (3 * real_mults_per_fft) + (N * 4)
    total_real_adds_fft = (3 * real_adds_per_fft) + (N * 2)

    time_fft_mult = total_real_mults_fft * float_mult_time
    time_fft_add = total_real_adds_fft * float_add_time
    total_time_fft = time_fft_mult + time_fft_add
    
    print(f"FFT size (N): {N}")
    print(f"Total floating point multiplications: {total_real_mults_fft}")
    print(f"Total floating point additions: {total_real_adds_fft}")
    print(f"Calculation: ({total_real_mults_fft} * {float_mult_time} ns) + ({total_real_adds_fft} * {float_add_time} ns)")
    print(f"Total time for FFT-based Convolution: {total_time_fft} ns\n")
    
    # --- Conclusion ---
    print("--- Comparison ---")
    times = {
        "Direct convolution with integers": total_time_direct_int,
        "Direct convolution with floating points": total_time_direct_float,
        "FFT-based convolution": total_time_fft,
    }
    
    fastest_method = min(times, key=times.get)
    fastest_time = times[fastest_method]

    print(f"Direct Convolution (Integers) time: {total_time_direct_int:,} ns")
    print(f"Direct Convolution (Floats) time:   {total_time_direct_float:,} ns")
    print(f"FFT-based Convolution time:       {total_time_fft:,} ns")
    print(f"\nThe fastest algorithm is '{fastest_method}' with an estimated time of {fastest_time:,} ns.")

if __name__ == "__main__":
    calculate_convolution_times()