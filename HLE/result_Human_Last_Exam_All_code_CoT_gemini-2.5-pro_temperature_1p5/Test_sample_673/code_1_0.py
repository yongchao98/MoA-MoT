import math

def solve():
    """
    Calculates and compares the execution time of two convolution algorithms.
    """
    # Given parameters
    n = 1000
    t_fp = 5  # ns per floating point operation
    t_int = 1  # ns per integer operation
    t_call = 15  # ns per function call

    # --- Algorithm 1: FFT-based ---
    print("Evaluating FFT-based Algorithm...")

    fft_num_calls = 2 * n - 1
    fft_fp_ops = 4 * n
    
    fft_time_from_calls = fft_num_calls * t_call
    fft_time_from_fp_ops = fft_fp_ops * t_fp
    fft_total_time = fft_time_from_calls + fft_time_from_fp_ops

    print(f"Time = (Number of function calls * Time per call) + (Number of FP ops * Time per FP op)")
    print(f"Time = ((2*n - 1) * {t_call}) + (4*n * {t_fp})")
    print(f"Time = ((2*{n} - 1) * {t_call}) + (4*{n} * {t_fp})")
    print(f"Time = ({fft_num_calls} * {t_call}) + ({fft_fp_ops} * {t_fp})")
    print(f"Time = {fft_time_from_calls} ns + {fft_time_from_fp_ops} ns")
    print(f"Total Time for FFT-based Algorithm = {fft_total_time} ns")
    print("-" * 30)

    # --- Algorithm 2: Direct with Fixed-Point ---
    print("Evaluating Direct Algorithm with Fixed-Point Arithmetic...")
    
    direct_fp_ops = 2 * n
    direct_int_ops = 2 * n**2
    # Assume the direct calculation is wrapped in a single function call
    direct_num_calls = 1

    direct_time_from_fp_ops = direct_fp_ops * t_fp
    direct_time_from_int_ops = direct_int_ops * t_int
    direct_time_from_calls = direct_num_calls * t_call
    direct_total_time = direct_time_from_fp_ops + direct_time_from_int_ops + direct_time_from_calls
    
    print(f"Time = (Number of FP ops * Time per FP op) + (Number of Int ops * Time per Int op) + (Number of function calls * Time per call)")
    print(f"Time = ((2*n) * {t_fp}) + ((2*n^2) * {t_int}) + (1 * {t_call})")
    print(f"Time = ((2*{n}) * {t_fp}) + ((2*{n}**2) * {t_int}) + (1 * {t_call})")
    print(f"Time = ({direct_fp_ops} * {t_fp}) + ({direct_int_ops} * {t_int}) + ({direct_num_calls} * {t_call})")
    print(f"Time = {direct_time_from_fp_ops} ns + {direct_time_from_int_ops} ns + {direct_time_from_calls} ns")
    print(f"Total Time for Direct Algorithm = {direct_total_time} ns")
    print("-" * 30)

    # --- Conclusion ---
    print("Comparison:")
    print(f"FFT-based Algorithm: {fft_total_time} ns")
    print(f"Direct Algorithm:      {direct_total_time} ns")

    if fft_total_time < direct_total_time:
        print("The original FFT-based algorithm is faster.")
        answer = "Y"
    else:
        print("The proposed direct algorithm is faster.")
        answer = "N"
    
    print("\nIs the original FFT-based algorithm is already the fastest?")
    print(f"<<<{answer}>>>")

solve()