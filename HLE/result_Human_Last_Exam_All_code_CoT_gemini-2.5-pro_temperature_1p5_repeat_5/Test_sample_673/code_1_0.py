import math

def solve():
    """
    Calculates and compares the execution times of two convolution algorithms.
    """
    # Given parameters
    n = 1000  # vector size
    t_float = 5  # ns per floating point operation
    t_int = 1  # ns per integer operation
    t_call = 15  # ns per function call

    # --- Algorithm 1: FFT-based ---
    # The 'divide-and-conquer' step implies recursive function calls. For a
    # size n problem, a recursive FFT implementation makes approximately 2n-1 calls.
    # The number of floating point operations is given as 4n.
    fft_calls = 2 * n - 1
    fft_fp_ops = 4 * n
    
    # Calculate each part of the cost for the equation
    fft_call_cost = fft_calls * t_call
    fft_fp_op_cost = fft_fp_ops * t_float
    total_time_fft = fft_call_cost + fft_fp_op_cost

    # --- Algorithm 2: Direct Integer Convolution ---
    # This involves 2n floating-point operations for conversion and 2n^2
    # integer operations for the direct convolution.
    direct_fp_ops = 2 * n
    direct_int_ops = 2 * n**2
    
    # Calculate each part of the cost for the equation
    direct_fp_op_cost = direct_fp_ops * t_float
    direct_int_op_cost = direct_int_ops * t_int
    total_time_direct = direct_fp_op_cost + direct_int_op_cost

    # --- Print the results ---
    print("Comparing algorithm speeds for n =", n)
    print("-" * 40)
    
    # Print detailed breakdown for FFT-based algorithm
    print("FFT-based Algorithm Time:")
    print(f"Equation: (2*n - 1)*T_call + (4*n)*T_float")
    print(f"Calculation: ({fft_calls}) * {t_call} ns + ({fft_fp_ops}) * {t_float} ns")
    print(f"Breakdown: {fft_call_cost} ns (from calls) + {fft_fp_op_cost} ns (from ops)")
    print(f"Total FFT Time = {total_time_fft} ns\n")

    # Print detailed breakdown for direct algorithm
    print("Direct Integer Algorithm Time:")
    print(f"Equation: (2*n)*T_float + (2*n^2)*T_int")
    print(f"Calculation: ({direct_fp_ops}) * {t_float} ns + ({direct_int_ops}) * {t_int} ns")
    print(f"Breakdown: {direct_fp_op_cost} ns (from conversion) + {direct_int_op_cost} ns (from convolution)")
    print(f"Total Direct Time = {total_time_direct} ns\n")
    
    # Final conclusion
    print("--- Conclusion ---")
    if total_time_fft < total_time_direct:
        print(f"The FFT-based algorithm ({total_time_fft} ns) is faster than the direct algorithm ({total_time_direct} ns).")
        # global answer
        # answer = "Y"
    else:
        print(f"The direct algorithm ({total_time_direct} ns) is faster than the FFT-based algorithm ({total_time_fft} ns).")
        # global answer
        # answer = "N"

solve()