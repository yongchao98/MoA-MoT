import math

def solve():
    """
    Calculates and compares the execution time of two convolution algorithms.
    """
    # 1. Define Constants from the problem description
    n = 1000
    t_float_op = 5  # Time in nanoseconds for a floating-point operation
    t_int_op = 1    # Time in nanoseconds for an integer operation
    t_func_call = 15 # Time in nanoseconds for a function call

    print("--- Algorithm 1: FFT-based ---")
    # 2. Calculate the total time for the FFT-based algorithm
    # We interpret the complexity as 4n*log2(n) float operations and 2n function calls.
    
    # Calculate operations and their time cost
    fft_float_ops = 4 * n * math.log2(n)
    fft_ops_time = fft_float_ops * t_float_op
    
    # Calculate function calls and their time cost
    fft_calls = 2 * n
    fft_calls_time = fft_calls * t_func_call
    
    # Calculate total time
    total_time_fft = fft_ops_time + fft_calls_time

    # Print the final equation with all numbers
    print("Equation: (Number of float ops * Time per op) + (Number of calls * Time per call)")
    print(f"Time = (4 * {n} * log2({n})) * {t_float_op} ns + (2 * {n}) * {t_func_call} ns")
    print(f"Time = ({fft_float_ops:.0f} ops * {t_float_op} ns) + ({fft_calls} calls * {t_func_call} ns)")
    print(f"Time = {fft_ops_time:.0f} ns + {fft_calls_time:.0f} ns = {total_time_fft:.0f} ns\n")


    print("--- Algorithm 2: Direct Integer-based ---")
    # 3. Calculate the total time for the direct integer-based algorithm
    
    # Calculate float operations for conversion and their time cost
    direct_float_ops = 2 * n
    direct_float_ops_time = direct_float_ops * t_float_op
    
    # Calculate integer operations for convolution and their time cost
    direct_int_ops = 2 * n**2
    direct_int_ops_time = direct_int_ops * t_int_op
    
    # Calculate total time
    total_time_direct = direct_float_ops_time + direct_int_ops_time

    # Print the final equation with all numbers
    print("Equation: (Num float ops * Time per op) + (Num int ops * Time per op)")
    print(f"Time = (2 * {n}) * {t_float_op} ns + (2 * {n}^2) * {t_int_op} ns")
    print(f"Time = ({direct_float_ops} ops * {t_float_op} ns) + ({direct_int_ops} ops * {t_int_op} ns)")
    print(f"Time = {direct_float_ops_time:.0f} ns + {direct_int_ops_time:.0f} ns = {total_time_direct:.0f} ns\n")
    

    print("--- Conclusion ---")
    # 4. Compare the two total times
    print(f"FFT-based algorithm time: \t{total_time_fft:.0f} ns")
    print(f"Direct integer algorithm time: \t{total_time_direct:.0f} ns")

    if total_time_fft < total_time_direct:
        print("\nThe FFT-based algorithm is faster.")
        final_answer = "Y"
    else:
        print("\nThe direct integer-based algorithm is faster.")
        final_answer = "N"
    
    print("\nIs the original FFT-based algorithm already the fastest? (Y/N)")
    print(f"<<<{final_answer}>>>")

solve()