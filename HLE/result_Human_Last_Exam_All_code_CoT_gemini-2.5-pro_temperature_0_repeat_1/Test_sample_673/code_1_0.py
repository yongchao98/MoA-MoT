import math

def solve():
    """
    This script compares the execution time of two convolution algorithms
    to determine which is faster for the given parameters.
    """

    # --- Given Parameters ---
    n = 1000  # Vector size
    t_float_op = 5  # Time for a floating-point operation in ns
    t_int_op = 1    # Time for an integer operation in ns
    t_func_call = 15 # Time for a function call in ns

    print("Step 1: Calculate the execution time for the direct convolution algorithm.\n")
    
    # --- Algorithm 2: Direct Convolution with Integers ---
    # This algorithm involves 2n^2 integer operations for the convolution,
    # 2n floating-point operations for conversions, and one function call.
    
    conv_ops_time = 2 * (n**2) * t_int_op
    conversion_time = 2 * n * t_float_op
    
    t_direct = conv_ops_time + conversion_time + t_func_call

    print("The total time for the direct algorithm is the sum of its parts:")
    print("Time_direct = (2 * n^2 * t_int_op) + (2 * n * t_float_op) + t_func_call")
    print(f"Time_direct = (2 * {n}^2 * {t_int_op}) + (2 * {n} * {t_float_op}) + {t_func_call}")
    print(f"Time_direct = {conv_ops_time} + {conversion_time} + {t_func_call}")
    print(f"Time_direct = {t_direct} ns\n")

    print("-" * 40)
    print("\nStep 2: Calculate the execution time for the FFT-based algorithm.\n")

    # --- Algorithm 1: FFT-based Algorithm ---
    # The time for this divide-and-conquer algorithm is modeled by the recurrence:
    # T(n) = 2 * T(n/2) + (cost at step n)
    # Cost at step n = 4n floating-point operations + 1 function call
    # The solved recurrence relation is approximately:
    # T(n) = (4 * t_float_op * log2(n) + 2 * t_func_call) * n - t_func_call
    
    log2_n = math.log2(n)
    
    # Using the derived closed-form solution from the recurrence relation
    # T(n) = (4*t_float*log2(n) + 2*t_call)*n - t_call
    # This is a good approximation for large n.
    term1_coeff = 4 * t_float_op
    term2_coeff = 2 * t_func_call
    
    t_fft = (term1_coeff * log2_n + term2_coeff) * n - t_func_call
    
    print("The time for the FFT-based algorithm is modeled by a recurrence relation.")
    print("An approximate closed-form solution for T(n) is:")
    print("Time_FFT = (4 * t_float_op * log2(n) + 2 * t_func_call) * n - t_func_call")
    print(f"Time_FFT = (4 * {t_float_op} * log2({n}) + 2 * {t_func_call}) * {n} - {t_func_call}")
    print(f"Time_FFT = ({term1_coeff} * {log2_n:.4f} + {term2_coeff}) * {n} - {t_func_call}")
    print(f"Time_FFT = ({term1_coeff * log2_n:.4f} + {term2_coeff}) * {n} - {t_func_call}")
    print(f"Time_FFT = {t_fft:.0f} ns\n")

    print("-" * 40)
    print("\nStep 3: Compare the two execution times.\n")
    print(f"FFT-based algorithm time: {t_fft:.0f} ns")
    print(f"Direct convolution time:  {t_direct} ns")

    if t_fft < t_direct:
        print("\nThe FFT-based algorithm is faster.")
        answer = "Y"
    else:
        print("\nThe direct convolution algorithm is faster.")
        answer = "N"
        
    print(f"\nIs the original FFT-based algorithm already the fastest? {answer}")
    print(f"\n<<<Y>>>")

solve()