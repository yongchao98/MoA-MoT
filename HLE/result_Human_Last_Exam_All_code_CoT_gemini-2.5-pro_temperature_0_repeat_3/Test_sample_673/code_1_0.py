import math

def solve():
    """
    Calculates and compares the execution times of two convolution algorithms.
    """
    # Step 1: Define the given parameters
    n = 1000  # Typical vector size
    t_float_op = 5  # ns per floating point operation
    t_int_op = 1    # ns per integer operation
    t_call = 15     # ns per function call

    # Step 2: Calculate the total time for the FFT-based algorithm
    print("--- FFT-based Algorithm Time Calculation ---")
    
    # Number of function calls for 3 FFTs (2 forward, 1 inverse)
    fft_calls = 3 * (2 * n - 1)
    fft_call_time = fft_calls * t_call
    
    # Number of floating point operations
    fft_float_ops = 4 * n
    fft_op_time = fft_float_ops * t_float_op
    
    # Total time for FFT-based algorithm
    total_time_fft = fft_call_time + fft_op_time
    
    print(f"Equation: T_FFT = (3 * (2 * n - 1)) * t_call + (4 * n) * t_float_op")
    print(f"Substituting values: T_FFT = (3 * (2 * {n} - 1)) * {t_call} + (4 * {n}) * {t_float_op}")
    print(f"Time for function calls = {fft_calls} calls * {t_call} ns/call = {fft_call_time} ns")
    print(f"Time for float operations = {fft_float_ops} ops * {t_float_op} ns/op = {fft_op_time} ns")
    print(f"Total T_FFT = {fft_call_time} ns + {fft_op_time} ns = {total_time_fft} ns\n")

    # Step 3: Calculate the total time for the direct integer-based algorithm
    print("--- Direct Integer-based Algorithm Time Calculation ---")
    
    # Time for floating point operations (conversion)
    direct_float_ops = 2 * n
    direct_conversion_time = direct_float_ops * t_float_op
    
    # Time for integer operations (convolution)
    direct_int_ops = 2 * n**2
    direct_convolution_time = direct_int_ops * t_int_op
    
    # Time for a single function call to the convolution routine
    direct_call_time = 1 * t_call
    
    # Total time for direct algorithm
    total_time_direct = direct_conversion_time + direct_convolution_time + direct_call_time
    
    print(f"Equation: T_Direct = (2 * n^2) * t_int_op + (2 * n) * t_float_op + 1 * t_call")
    print(f"Substituting values: T_Direct = (2 * {n}^2) * {t_int_op} + (2 * {n}) * {t_float_op} + 1 * {t_call}")
    print(f"Time for integer convolution = {direct_int_ops} ops * {t_int_op} ns/op = {direct_convolution_time} ns")
    print(f"Time for float conversion = {direct_float_ops} ops * {t_float_op} ns/op = {direct_conversion_time} ns")
    print(f"Time for function call = 1 call * {t_call} ns/call = {direct_call_time} ns")
    print(f"Total T_Direct = {direct_convolution_time} ns + {direct_conversion_time} ns + {direct_call_time} ns = {total_time_direct} ns\n")
    
    # Step 4: Compare the two total times
    print("--- Comparison ---")
    print(f"FFT-based algorithm time: {total_time_fft} ns")
    print(f"Direct integer algorithm time: {total_time_direct} ns")
    
    is_fft_faster = total_time_fft < total_time_direct
    answer = 'Y' if is_fft_faster else 'N'
    
    print(f"\nIs the original FFT-based algorithm already the fastest? ({total_time_fft} < {total_time_direct})")
    print(f"{answer}")
    
    print(f"<<<{answer}>>>")

solve()