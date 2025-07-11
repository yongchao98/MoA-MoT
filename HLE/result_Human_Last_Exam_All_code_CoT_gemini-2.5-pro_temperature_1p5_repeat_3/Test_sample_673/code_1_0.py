import math

def solve():
    """
    Compares the execution time of two convolution algorithms.
    """
    # Given parameters
    n = 1000
    T_float = 5  # ns per floating point operation
    T_int = 1    # ns per integer operation
    T_call = 15  # ns per function call

    # --- Algorithm 1: FFT-based ---
    # Cost of floating point operations
    fft_fp_ops = 4 * n
    fft_fp_time = fft_fp_ops * T_float
    
    # Cost of divide-and-conquer step (modeled as function call overhead)
    # A recursive FFT of size n makes 2n-1 calls
    fft_calls = 2 * n - 1
    fft_call_time = fft_calls * T_call
    
    # Total time for FFT-based algorithm
    T_fft_total = fft_fp_time + fft_call_time

    print("--- FFT-based Algorithm ---")
    print(f"Time for floating point operations = {fft_fp_ops} ops * {T_float} ns/op = {fft_fp_time} ns")
    print(f"Time for function call overhead = {fft_calls} calls * {T_call} ns/call = {fft_call_time} ns")
    print(f"Total FFT Time = {fft_fp_time} ns + {fft_call_time} ns = {T_fft_total} ns")
    print("\nIn one equation:")
    print(f"T_FFT = (4 * {n} * {T_float}) + ((2 * {n} - 1) * {T_call}) = {T_fft_total} ns\n")


    # --- Algorithm 2: Direct Convolution ---
    # Cost of float-to-int conversion
    direct_fp_ops = 2 * n
    direct_fp_time = direct_fp_ops * T_float
    
    # Cost of integer convolution
    direct_int_ops = 2 * n**2
    direct_int_time = direct_int_ops * T_int

    # Total time for direct convolution algorithm
    T_direct_total = direct_fp_time + direct_int_time
    
    print("--- Direct Convolution Algorithm ---")
    print(f"Time for data conversion = {direct_fp_ops} ops * {T_float} ns/op = {direct_fp_time} ns")
    print(f"Time for integer convolution = {direct_int_ops} ops * {T_int} ns/op = {direct_int_time} ns")
    print(f"Total Direct Convolution Time = {direct_fp_time} ns + {direct_int_time} ns = {T_direct_total} ns")
    print("\nIn one equation:")
    print(f"T_Direct = (2 * {n} * {T_float}) + (2 * {n}^2 * {T_int}) = {T_direct_total} ns\n")


    # --- Comparison ---
    print("--- Conclusion ---")
    if T_fft_total < T_direct_total:
        print(f"{T_fft_total} ns (FFT) < {T_direct_total} ns (Direct)")
        print("The original FFT-based algorithm is already the fastest.")
        answer = "Y"
    else:
        print(f"{T_fft_total} ns (FFT) >= {T_direct_total} ns (Direct)")
        print("The proposed direct convolution algorithm is faster.")
        answer = "N"

    print(f"\nIs the original FFT-based algorithm is already the fastest? \nAnswer: {answer}")
    print(f"<<<{answer}>>>")

solve()