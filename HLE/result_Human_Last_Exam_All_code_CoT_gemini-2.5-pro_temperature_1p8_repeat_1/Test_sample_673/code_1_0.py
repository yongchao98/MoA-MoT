import math

def solve_algorithm_comparison():
    """
    This function calculates and compares the execution times for two convolution algorithms
    based on the parameters given in the problem description.
    """
    # Given parameters
    n = 1000
    t_float_op = 5  # ns
    t_int_op = 1    # ns
    t_call = 15     # ns

    print("Analyzing the performance of two convolution algorithms for vector size n = 1000.")
    print(f"Unit costs: Floating point op = {t_float_op} ns, Integer op = {t_int_op} ns, Function call = {t_call} ns\n")

    # --- Algorithm 1: FFT-based ---
    print("--- Algorithm 1: FFT-based Calculation ---")

    # Cost from the divide-and-conquer step (function calls)
    num_calls_fft = 2 * n - 1
    cost_calls_fft = num_calls_fft * t_call

    # Cost from floating point operations
    num_float_ops_fft = 4 * n
    cost_float_ops_fft = num_float_ops_fft * t_float_op

    # Total time for FFT-based algorithm
    total_time_fft = cost_calls_fft + cost_float_ops_fft

    print(f"Equation for the total time: ((2 * n - 1) * t_call) + (4 * n * t_float)")
    print(f"Time for divide-and-conquer step: ({2} * {n} - 1) * {t_call} = {num_calls_fft} * {t_call} = {cost_calls_fft} ns")
    print(f"Time for floating point operations: {4} * {n} * {t_float_op} = {num_float_ops_fft} * {t_float_op} = {cost_float_ops_fft} ns")
    print(f"Total time = {cost_calls_fft} + {cost_float_ops_fft} = {total_time_fft} ns\n")


    # --- Algorithm 2: Direct convolution with fixed-point integers ---
    print("--- Algorithm 2: Direct Integer Convolution Calculation ---")

    # Cost from conversion to fixed-point
    num_float_ops_direct = 2 * n
    cost_conversion_direct = num_float_ops_direct * t_float_op

    # Cost from integer convolution
    num_int_ops_direct = 2 * n**2
    cost_convolution_direct = num_int_ops_direct * t_int_op

    # Total time for direct convolution algorithm
    total_time_direct = cost_conversion_direct + cost_convolution_direct
    
    print(f"Equation for the total time: (2 * n * t_float) + (2 * n^2 * t_int)")
    print(f"Time for conversion: {2} * {n} * {t_float_op} = {num_float_ops_direct} * {t_float_op} = {cost_conversion_direct} ns")
    print(f"Time for integer convolution: {2} * {n}^2 * {t_int_op} = {f'{num_int_ops_direct:,}'} * {t_int_op} = {f'{cost_convolution_direct:,}'} ns")
    print(f"Total time = {cost_conversion_direct} + {cost_convolution_direct} = {f'{total_time_direct:,}'} ns\n")

    # --- Comparison ---
    print("--- Comparison ---")
    print(f"FFT-based algorithm time: {f'{total_time_fft:,}'} ns")
    print(f"Direct convolution time: {f'{total_time_direct:,}'} ns")

    if total_time_fft < total_time_direct:
        print("\nThe original FFT-based algorithm is faster.")
        answer = "Y"
    else:
        print("\nThe direct convolution algorithm is faster.")
        answer = "N"

    print(f"\n<<<{answer}>>>")

# Execute the function to solve the problem
solve_algorithm_comparison()