import math

def solve():
    """
    Calculates and compares the execution time of two convolution algorithms.
    """
    # --- Given Parameters ---
    n = 1000
    t_float = 5  # ns per floating point operation
    t_int = 1    # ns per integer operation
    t_call = 15  # ns per function call

    print("--- Analysis of Algorithm 1: FFT-based ---")
    print("This algorithm is recursive (divide-and-conquer).")
    print("Time = (Function Call Cost) + (Operations Cost)")
    
    # --- Time for FFT-based Algorithm (T1) ---
    # The time is modeled by the equation: T_fft = (2*n - 1) * t_call + 4 * n * log2(n) * t_float
    print(f"\nEquation: T_fft = (2 * n - 1) * t_call + 4 * n * log2(n) * t_float")

    # Cost of function calls
    num_calls_fft = 2 * n - 1
    time_calls_fft = num_calls_fft * t_call

    # Cost of floating point operations
    log2_n = math.log2(n)
    num_ops_fft_total = 4 * n * log2_n
    time_ops_fft = num_ops_fft_total * t_float

    # Total time for FFT-based algorithm
    total_time_fft = time_calls_fft + time_ops_fft

    print(f"\nCalculation for n = {n}:")
    print(f"  Function Call Cost = ({2} * {n} - 1) * {t_call} ns = {num_calls_fft} * {t_call} ns = {time_calls_fft} ns")
    print(f"  Operations Cost = {4} * {n} * log2({n}) * {t_float} ns = {4 * n} * {log2_n:.4f} * {t_float} ns = {time_ops_fft:.2f} ns")
    print(f"Total time for FFT-based algorithm (T1) = {time_calls_fft} ns + {time_ops_fft:.2f} ns = {total_time_fft:.2f} ns")

    print("\n" + "="*50 + "\n")

    print("--- Analysis of Algorithm 2: Integer-based Direct Convolution ---")
    print("This algorithm is direct (non-recursive).")
    print("Time = (Function Call Cost) + (Conversion Cost) + (Convolution Cost)")
    
    # --- Time for Integer-based Algorithm (T2) ---
    # The time is modeled by the equation: T_direct = 1 * t_call + 2 * n * t_float + 2 * n^2 * t_int
    print(f"\nEquation: T_direct = 1 * t_call + 2 * n * t_float + 2 * n^2 * t_int")

    # Cost of the single function call
    time_call_direct = 1 * t_call

    # Cost of float-to-int conversion
    num_ops_conv = 2 * n
    time_ops_conv = num_ops_conv * t_float

    # Cost of integer convolution
    num_ops_conv_int = 2 * n**2
    time_ops_conv_int = num_ops_conv_int * t_int

    # Total time for integer-based algorithm
    total_time_direct = time_call_direct + time_ops_conv + time_ops_conv_int
    
    print(f"\nCalculation for n = {n}:")
    print(f"  Function Call Cost = {1} * {t_call} ns = {time_call_direct} ns")
    print(f"  Conversion Cost = {2} * {n} * {t_float} ns = {num_ops_conv} * {t_float} ns = {time_ops_conv} ns")
    print(f"  Convolution Cost = {2} * {n}^2 * {t_int} ns = {num_ops_conv_int} * {t_int} ns = {time_ops_conv_int} ns")
    print(f"Total time for integer-based algorithm (T2) = {time_call_direct} ns + {time_ops_conv} ns + {time_ops_conv_int} ns = {total_time_direct:.2f} ns")

    print("\n" + "="*50 + "\n")

    print("--- Comparison and Conclusion ---")
    print(f"T1 (FFT-based) = {total_time_fft:.2f} ns")
    print(f"T2 (Integer-based) = {total_time_direct:.2f} ns")

    if total_time_fft < total_time_direct:
        print("\nThe original FFT-based algorithm is faster.")
        answer = "Y"
    else:
        print("\nThe integer-based direct convolution algorithm is faster.")
        answer = "N"
    
    print(f"\nIs the original FFT-based algorithm is already the fastest?")
    print(f"<<<{answer}>>>")

solve()