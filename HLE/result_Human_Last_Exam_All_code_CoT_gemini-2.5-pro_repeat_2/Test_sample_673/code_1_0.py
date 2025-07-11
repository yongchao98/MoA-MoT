import math

def solve():
    """
    Calculates and compares the execution time of two convolution algorithms
    to determine which one is faster.
    """
    # --- Step 1: Define Constants ---
    n = 1000.0  # vector size
    t_float = 5   # ns per floating-point operation
    t_int = 1     # ns per integer operation
    t_call = 15   # ns per function call

    print("To determine the fastest algorithm, we calculate the total execution time for each.")
    print(f"Parameters: n={int(n)}, t_float={t_float}ns, t_int={t_int}ns, t_call={t_call}ns\n")

    # --- Step 2: Calculate Time for Algorithm 1 (FFT-based) ---
    print("--- Algorithm 1: FFT-based ---")
    log2_n = math.log2(n)

    # Cost of the "divide-and-conquer step" (3 FFTs)
    fft_ops = 3 * (5 * n * log2_n)
    fft_calls = 3 * (2 * n)
    time_d_and_c = (fft_ops * t_float) + (fft_calls * t_call)

    # Cost of the "final calculation"
    final_ops = 4 * n
    time_final_calc = final_ops * t_float

    # Total time for Algorithm 1
    time_algo1 = time_d_and_c + time_final_calc

    print("Equation: (3*(5*n*log2(n)) * t_float + 3*(2*n) * t_call) + (4*n * t_float)")
    print(f"Time = (3*(5*{int(n)}*log2({int(n)})) * {t_float} + 3*(2*{int(n)}) * {t_call}) + (4*{int(n)} * {t_float})")
    print(f"Time = (3*(5*{int(n)}*{log2_n:.2f}) * {t_float} + {int(fft_calls)} * {t_call}) + ({int(final_ops)} * {t_float})")
    print(f"Time = ({fft_ops:.0f} * {t_float} + {int(fft_calls)} * {t_call}) + {int(time_final_calc)}")
    print(f"Time = {time_d_and_c:.0f} + {time_final_calc:.0f}")
    print(f"Total Time for FFT-based algorithm = {time_algo1:.0f} ns\n")


    # --- Step 3: Calculate Time for Algorithm 2 (Integer-based) ---
    print("--- Algorithm 2: Integer-based ---")
    conversion_ops = 2 * n
    convolution_ops = 2 * n**2
    
    time_conversion = conversion_ops * t_float
    time_convolution = convolution_ops * t_int
    time_algo2 = time_conversion + time_convolution

    print("Equation: (2*n * t_float) + (2*n^2 * t_int)")
    print(f"Time = (2*{int(n)} * {t_float}) + (2*{int(n)}^2 * {t_int})")
    print(f"Time = ({int(time_conversion)}) + ({int(time_convolution)})")
    print(f"Total Time for Integer-based algorithm = {time_algo2:.0f} ns\n")

    # --- Step 4: Compare and Conclude ---
    print("--- Conclusion ---")
    if time_algo1 < time_algo2:
        print(f"{time_algo1:.0f} ns (FFT) < {time_algo2:.0f} ns (Integer). The original algorithm is faster.")
        answer = 'Y'
    else:
        print(f"{time_algo1:.0f} ns (FFT) >= {time_algo2:.0f} ns (Integer). The original algorithm is not the fastest.")
        answer = 'N'

    print("\nIs the original FFT-based algorithm already the fastest?")
    print(f"<<<{answer}>>>")

solve()