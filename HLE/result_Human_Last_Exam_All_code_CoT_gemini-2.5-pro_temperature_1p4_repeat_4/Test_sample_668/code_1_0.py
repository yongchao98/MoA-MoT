import math

def calculate_convolution_times():
    """
    Calculates and compares the time for different convolution algorithms
    on a specific machine.
    """
    # Problem parameters
    n = 1000
    t_int_add = 1  # ns
    t_int_mul = 2  # ns
    t_fp_add = 9   # ns
    t_fp_mul = 19  # ns

    print("Determining the fastest convolution algorithm for n=1000.\n")
    print(f"Machine Operation Times:\n- Integer Add/Sub: {t_int_add} ns\n- Integer Mul: {t_int_mul} ns\n- Floating Point Add/Sub: {t_fp_add} ns\n- Floating Point Mul: {t_fp_mul} ns\n")
    print("-" * 60)

    # --- Method B: Direct convolution with integers ---
    print("Method B: Direct convolution with integers")
    # Number of operations for direct convolution of two n-element series
    num_mul_direct = n * n
    num_add_direct = n * n - n

    # Time calculation
    time_mul_int = num_mul_direct * t_int_mul
    time_add_int = num_add_direct * t_int_add
    total_time_direct_int = time_mul_int + time_add_int

    print(f"The number of multiplications is n^2 = {n} * {n} = {num_mul_direct}")
    print(f"The number of additions is n^2 - n = {n}*{n} - {n} = {num_add_direct}")
    print("Equation for total time:")
    print(f"Total time = (Number of multiplications * Integer multiplication time) + (Number of additions * Integer addition time)")
    print(f"Total time = ({num_mul_direct} * {t_int_mul} ns) + ({num_add_direct} * {t_int_add} ns)")
    print(f"Total time = {time_mul_int} ns + {time_add_int} ns = {total_time_direct_int} ns")
    print("-" * 60)

    # --- Method C: Direct convolution with floating points ---
    print("Method C: Direct convolution with floating points")
    # Number of operations is the same as for integers
    time_mul_fp = num_mul_direct * t_fp_mul
    time_add_fp = num_add_direct * t_fp_add
    total_time_direct_fp = time_mul_fp + time_add_fp

    print(f"The number of multiplications is n^2 = {n} * {n} = {num_mul_direct}")
    print(f"The number of additions is n^2 - n = {n}*{n} - {n} = {num_add_direct}")
    print("Equation for total time:")
    print(f"Total time = (Number of multiplications * Floating point multiplication time) + (Number of additions * Floating point addition time)")
    print(f"Total time = ({num_mul_direct} * {t_fp_mul} ns) + ({num_add_direct} * {t_fp_add} ns)")
    print(f"Total time = {time_mul_fp} ns + {time_add_fp} ns = {total_time_direct_fp} ns")
    print("-" * 60)

    # --- Method A: FFT-based Convolution ---
    print("Method A: FFT-based convolution")
    # Step 1: Find N, the padded size. It must be a power of 2 and >= 2n-1.
    padded_len_req = 2 * n - 1
    N_power = math.ceil(math.log2(padded_len_req))
    N = int(2**N_power)
    log2_N = N_power

    print(f"The length of the convolved sequence is 2*n-1 = 2*{n}-1 = {padded_len_req}.")
    print(f"Data is padded to the next power of 2, so the FFT size N = {N}.")
    print(f"log2(N) = log2({N}) = {int(log2_N)}.")

    # Number of operations for the entire process (2x FFT, 1x complex mult, 1x IFFT)
    # A standard complex N-point FFT requires ~2*N*log2(N) real multiplications and ~3*N*log2(N) real additions.
    # Total real multiplications: 2xFFT + N complex mult + 1xIFFT + N scaling
    num_mul_fft_ffts = 2 * (2 * N * log2_N)
    num_mul_fft_elementwise = N * 4
    num_mul_fft_ifft = 1 * (2 * N * log2_N)
    num_mul_fft_scaling = N
    num_mul_fft_total = num_mul_fft_ffts + num_mul_fft_elementwise + num_mul_fft_ifft + num_mul_fft_scaling

    # Total real additions: 2xFFT + N complex mult + 1xIFFT
    num_add_fft_ffts = 2 * (3 * N * log2_N)
    num_add_fft_elementwise = N * 2
    num_add_fft_ifft = 1 * (3 * N * log2_N)
    num_add_fft_total = num_add_fft_ffts + num_add_fft_elementwise + num_add_fft_ifft
    
    print("Floating point operations breakdown:")
    print(f"Multiplications = (2*FFTs) + (Element-wise prod) + (IFFT) + (IFFT scaling)")
    print(f"Number of mults = {int(num_mul_fft_ffts)} + {int(num_mul_fft_elementwise)} + {int(num_mul_fft_ifft)} + {int(num_mul_fft_scaling)} = {int(num_mul_fft_total)}")
    print(f"Additions = (2*FFTs) + (Element-wise prod) + (IFFT)")
    print(f"Number of adds = {int(num_add_fft_ffts)} + {int(num_add_fft_elementwise)} + {int(num_add_fft_ifft)} = {int(num_add_fft_total)}")

    # Time calculation for FFT
    time_mul_fft = num_mul_fft_total * t_fp_mul
    time_add_fft = num_add_fft_total * t_fp_add
    total_time_fft = time_mul_fft + time_add_fft

    print("Equation for total time:")
    print(f"Total time = ({int(num_mul_fft_total)} * {t_fp_mul} ns) + ({int(num_add_fft_total)} * {t_fp_add} ns)")
    print(f"Total time = {int(time_mul_fft)} ns + {int(time_add_fft)} ns = {int(total_time_fft)} ns")
    print("-" * 60)

    # --- Conclusion ---
    print("Comparison and Conclusion")
    print(f"Time for Method A (FFT):                       {int(total_time_fft)} ns")
    print(f"Time for Method B (Direct, Integers):          {total_time_direct_int} ns")
    print(f"Time for Method C (Direct, Floating Point):    {total_time_direct_fp} ns")
    print("-" * 60)

    times = {
        'A': total_time_fft,
        'B': total_time_direct_int,
        'C': total_time_direct_fp,
    }

    fastest_method_key = min(times, key=times.get)
    fastest_method_name = {
        'A': 'FFT',
        'B': 'Direct convolution with integers',
        'C': 'Direct convolution with floating points'
    }[fastest_method_key]

    print(f"The fastest algorithm is '{fastest_method_name}' (Method {fastest_method_key}).")


if __name__ == '__main__':
    calculate_convolution_times()