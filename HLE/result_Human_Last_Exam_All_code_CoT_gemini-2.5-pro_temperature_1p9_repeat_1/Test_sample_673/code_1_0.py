import math

def solve_convolution_speed():
    """
    Calculates and compares the execution time for two convolution algorithms.
    """
    # --- Given Parameters ---
    n = 1000  # vector size
    t_float = 5  # ns per floating point operation
    t_int = 1    # ns per integer operation
    t_call = 15  # ns per function call

    # --- Print Preamble and Parameters ---
    print("Comparing two convolution algorithms for a vector of size n")
    print("-" * 60)
    print("Parameters:")
    print(f"  Vector size (n): {n}")
    print(f"  Floating point operation time (t_float): {t_float} ns")
    print(f"  Integer operation time (t_int): {t_int} ns")
    print(f"  Function call time (t_call): {t_call} ns")
    print("-" * 60)

    # --- Algorithm 1: FFT-based Method ---
    # Cost is modeled by the recurrence T(n) = 2*T(n/2) + (4n*t_float + t_call).
    # Solved, this gives T(n) = (4 * n * log2(n) * t_float) + ((2n-1) * t_call).
    log2_n = math.log2(n)
    
    # Calculate each term in the cost formula
    fft_work_term_val = 4 * n * log2_n * t_float
    fft_call_term_val = (2 * n - 1) * t_call
    time_fft = fft_work_term_val + fft_call_term_val

    print("Algorithm 1: FFT-based Method")
    print("Cost Formula: (4 * n * log2(n) * t_float) + ((2*n - 1) * t_call)")
    print(f"Calculation: (4 * {n} * log2({n}) * {t_float}) + ((2*{n} - 1) * {t_call})")
    print(f"             = (4 * {n} * {log2_n:.4f} * {t_float}) + ({2*n-1} * {t_call})")
    print(f"             = {fft_work_term_val:.2f} + {fft_call_term_val:.2f}")
    print(f"Total time (FFT): {time_fft:.2f} ns")
    print("-" * 60)

    # --- Algorithm 2: Direct Method with Fixed-Point Integers ---
    # Cost is the sum of its components: float conversion, integer convolution, and one function call.
    # Cost Formula: (2 * n^2 * t_int) + (2 * n * t_float) + t_call
    
    # Calculate each term in the cost formula
    direct_conv_term_val = 2 * n**2 * t_int
    direct_conversion_term_val = 2 * n * t_float
    direct_call_term_val = t_call
    time_direct = direct_conv_term_val + direct_conversion_term_val + direct_call_term_val

    print("Algorithm 2: Direct Method with Fixed-Point Integers")
    print("Cost Formula: (2 * n^2 * t_int) + (2 * n * t_float) + t_call")
    print(f"Calculation: (2 * {n}^2 * {t_int}) + (2 * {n} * {t_float}) + {t_call}")
    print(f"             = (2 * {n**2} * {t_int}) + (2 * {n} * {t_float}) + {t_call}")
    print(f"             = {direct_conv_term_val:.2f} + {direct_conversion_term_val:.2f} + {direct_call_term_val:.2f}")
    print(f"Total time (Direct): {time_direct:.2f} ns")
    print("-" * 60)
    
    # --- Conclusion ---
    print("Conclusion:")
    if time_fft < time_direct:
        print(f"The FFT-based algorithm ({time_fft:.0f} ns) is faster than the direct algorithm ({time_direct:.0f} ns).")
        answer = "Y"
    else:
        print(f"The direct algorithm ({time_direct:.0f} ns) is faster than the FFT-based algorithm ({time_fft:.0f} ns).")
        answer = "N"

    print("\nIs the original FFT-based algorithm is already the fastest? (Y/N)")
    print(f"<<<{answer}>>>")

solve_convolution_speed()