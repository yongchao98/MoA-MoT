import math

def solve_algorithm_speed():
    """
    Compares the speed of two convolution algorithms based on provided parameters.
    """
    # Given parameters
    t_fp_op = 5   # Time for a floating-point operation in ns
    t_int_op = 1  # Time for an integer operation in ns
    t_call = 15   # Time for a function call in ns
    n = 1000      # Typical vector size

    print("--- Analysis of Algorithm Speeds ---")
    print(f"Given parameters:\n- Vector size (n): {n}\n- Float op time: {t_fp_op} ns\n- Integer op time: {t_int_op} ns\n- Function call time: {t_call} ns\n")

    # --- Algorithm 1: FFT-based (Original) ---
    # The time is modeled by the recurrence T(n) = 2T(n/2) + cost_of_combination.
    # For a recursive FFT, this solves to (2n-1) function calls and 4n*log2(n) floating point operations.
    log2_n = math.log2(n)
    
    # Cost from recursive function calls
    fft_call_cost = (2 * n - 1) * t_call
    # Cost from floating point operations
    fft_op_cost = (4 * n * log2_n) * t_fp_op
    
    time_fft = fft_call_cost + fft_op_cost

    print("--- Algorithm 1: FFT-based ---")
    print("Time Formula = (Number of Calls * Time per Call) + (Number of FP Ops * Time per FP Op)")
    print(f"Equation: (2 * {n} - 1) * {t_call} + (4 * {n} * log2({n})) * {t_fp_op}")
    print(f"Number of Calls       = 2 * {n} - 1 = {2 * n - 1}")
    print(f"Number of FP Ops      = 4 * {n} * {log2_n:.4f} = {4 * n * log2_n:.2f}")
    print(f"Time (FFT) = ({2*n-1} * {t_call}) + ({4*n*log2_n:.2f} * {t_fp_op})")
    print(f"           = {fft_call_cost} ns (calls) + {fft_op_cost:.2f} ns (ops)")
    print(f"Total Time = {time_fft:.2f} ns\n")


    # --- Algorithm 2: Direct Convolution (Fixed-Point) ---
    # This involves 1 function call, 2n float->int conversions, and 2n^2 integer operations.
    
    # Cost from single function call
    direct_call_cost = t_call
    # Cost from converting 2 vectors from float to integer
    direct_conversion_cost = (2 * n) * t_fp_op
    # Cost from direct integer convolution
    direct_conv_cost = (2 * n**2) * t_int_op
    
    time_direct = direct_call_cost + direct_conversion_cost + direct_conv_cost

    print("--- Algorithm 2: Direct Convolution ---")
    print("Time Formula = Time per Call + (Conversion Ops * Time per FP Op) + (Conv Ops * Time per Int Op)")
    print(f"Equation: {t_call} + (2 * {n}) * {t_fp_op} + (2 * {n}^2) * {t_int_op}")
    print(f"Number of FP Ops (Conversion) = 2 * {n} = {2 * n}")
    print(f"Number of Int Ops (Convolution) = 2 * {n}^2 = {2 * n**2}")
    print(f"Time (Direct) = {t_call} + ({2 * n} * {t_fp_op}) + ({2 * n**2} * {t_int_op})")
    print(f"              = {direct_call_cost} ns (call) + {direct_conversion_cost} ns (conversion) + {direct_conv_cost} ns (convolution)")
    print(f"Total Time    = {time_direct:.2f} ns\n")


    # --- Conclusion ---
    print("--- Conclusion ---")
    if time_fft < time_direct:
        print(f"Comparing the total times ({time_fft:.0f} ns vs {time_direct:.0f} ns), the original FFT-based algorithm is faster.")
        answer = "Y"
    else:
        print(f"Comparing the total times ({time_fft:.0f} ns vs {time_direct:.0f} ns), the direct convolution algorithm is faster.")
        answer = "N"

    print(f"\nIs the original FFT-based algorithm already the fastest? {answer}")

solve_algorithm_speed()