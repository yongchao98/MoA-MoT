import math

def solve_convolution_time():
    """
    Calculates and compares the execution time for two convolution algorithms.
    """
    # Step 1 & 2: Define constants and parameters
    n = 1000
    t_float = 5  # ns per operation
    t_int = 1    # ns per operation
    t_call = 15  # ns per call

    print("Comparing two convolution algorithms for a vector of size n = 1000.\n")
    print("--- Given Parameters ---")
    print(f"Vector size (n): {n}")
    print(f"Floating point operation time (t_float): {t_float} ns")
    print(f"Integer operation time (t_int): {t_int} ns")
    print(f"Function call time (t_call): {t_call} ns\n")

    # --- Algorithm 1: FFT-based Method ---
    print("--- 1. FFT-based Algorithm ---")
    print("The time complexity is based on the recurrence T(n) = 2*T(n/2) + (4*n*t_float + t_call).")
    print("The solved closed-form equation for the time is: (2 * t_call * n) + (4 * t_float * n * log2(n)) - t_call\n")

    log_n_val = math.log2(n)
    # Perform the calculation
    term1_fft = 2 * t_call * n
    term2_fft = 4 * t_float * n * log_n_val
    term3_fft = -t_call
    time_fft = term1_fft + term2_fft + term3_fft

    print("Equation with values:")
    print(f"(2 * {t_call}) * {n} + (4 * {t_float}) * {n} * log2({n}) - {t_call}")
    print(f"= {term1_fft} + {term2_fft:.0f} - {abs(term3_fft)}")
    print(f"= {time_fft:.0f} ns\n")


    # --- Algorithm 2: Direct Integer-based Method ---
    print("--- 2. Direct Integer-based Algorithm ---")
    print("The time is the sum of conversion, convolution, and function call costs.")
    print("Equation for the time is: (2 * n * t_float) + (2 * n^2 * t_int) + t_call\n")

    # Perform the calculation
    term1_direct = 2 * n * t_float
    term2_direct = 2 * (n**2) * t_int
    term3_direct = t_call
    time_direct = term1_direct + term2_direct + term3_direct

    print("Equation with values:")
    print(f"(2 * {n} * {t_float}) + (2 * {n}^2 * {t_int}) + {t_call}")
    print(f"= {term1_direct} + {term2_direct} + {term3_direct}")
    print(f"= {time_direct:.0f} ns\n")


    # Step 4: Compare and conclude
    print("--- Conclusion ---")
    print(f"FFT-based algorithm time: {time_fft:.0f} ns")
    print(f"Direct algorithm time:      {time_direct:.0f} ns")

    if time_fft < time_direct:
        print("\nThe original FFT-based algorithm is faster.")
        answer = "Y"
    else:
        print("\nThe direct integer-based algorithm is faster.")
        answer = "N"

    print(f"<<<{answer}>>>")

# Run the analysis
solve_convolution_time()