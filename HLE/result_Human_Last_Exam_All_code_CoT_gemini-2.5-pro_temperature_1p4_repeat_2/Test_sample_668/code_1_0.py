import math

def solve():
    """
    Calculates the fastest convolution algorithm for a given machine and series size.
    """
    # Machine operation times in nanoseconds (ns)
    t_int_add = 1
    t_int_mul = 2
    t_float_add = 9
    t_float_mul = 19

    # Series size
    n = 1000

    print("Step 1: Time Estimation for Direct Convolution\n")
    
    # --- Direct convolution with Integers ---
    num_int_muls = n ** 2
    num_int_adds = (n - 1) ** 2
    time_int_muls = num_int_muls * t_int_mul
    time_int_adds = num_int_adds * t_int_add
    total_time_int = time_int_muls + time_int_adds
    
    print("Direct Convolution with Integers:")
    print(f"Number of multiplications: {n}^2 = {num_int_muls}")
    print(f"Number of additions: ({n}-1)^2 = {num_int_adds}")
    print(f"Total time = ({num_int_muls} * {t_int_mul}) + ({num_int_adds} * {t_int_add}) = {time_int_muls} + {time_int_adds} = {total_time_int} ns\n")

    # --- Direct convolution with Floating Points ---
    num_float_muls = n ** 2
    num_float_adds = (n - 1) ** 2
    time_float_muls = num_float_muls * t_float_mul
    time_float_adds = num_float_adds * t_float_add
    total_time_float = time_float_muls + time_float_adds
    
    print("Direct Convolution with Floating Points:")
    print(f"Number of multiplications: {n}^2 = {num_float_muls}")
    print(f"Number of additions: ({n}-1)^2 = {num_float_adds}")
    print(f"Total time = ({num_float_muls} * {t_float_mul}) + ({num_float_adds} * {t_float_add}) = {time_float_muls} + {time_float_adds} = {total_time_float} ns\n")

    print("Step 2: Time Estimation for FFT-based Convolution\n")
    
    # --- FFT-based Convolution ---
    # Find the next power of 2 for FFT size N, where N >= 2n - 1
    required_len = 2 * n - 1
    N = 1
    while N < required_len:
        N *= 2
    log2N = int(math.log2(N))

    # Number of floating point operations for FFT
    # Real multiplications = 3 * (2*N*log2N) from 3 FFTs + 4*N from complex multiplication
    num_fft_muls = (6 * N * log2N) + (4 * N)
    # Real additions = 3 * (3*N*log2N) from 3 FFTs + 2*N from complex multiplication
    num_fft_adds = (9 * N * log2N) + (2 * N)
    
    time_fft_muls = num_fft_muls * t_float_mul
    time_fft_adds = num_fft_adds * t_float_add
    total_time_fft = time_fft_muls + time_fft_adds

    print("FFT-based Convolution (Floating Point):")
    print(f"FFT size N (power of 2 >= {required_len}): {N}")
    print(f"Number of multiplications: (6 * {N} * {log2N}) + (4 * {N}) = {num_fft_muls}")
    print(f"Number of additions: (9 * {N} * {log2N}) + (2 * {N}) = {num_fft_adds}")
    print(f"Total time = ({num_fft_muls} * {t_float_mul}) + ({num_fft_adds} * {t_float_add}) = {time_fft_muls} + {time_fft_adds} = {total_time_fft} ns\n")
    
    print("Step 3: Comparison\n")

    # --- Comparison ---
    results = {
        "Direct convolution with integers": total_time_int,
        "FFT-based convolution": total_time_fft,
        "Direct convolution with floating points": total_time_float
    }
    
    fastest_method = min(results, key=results.get)
    
    print(f"Summary of total times:")
    for method, time in results.items():
        print(f"- {method}: {time} ns")
        
    print(f"\nThe fastest algorithm is: {fastest_method}")

if __name__ == '__main__':
    solve()