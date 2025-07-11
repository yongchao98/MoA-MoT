import math

def solve():
    """
    Calculates and compares the execution time for two convolution algorithms.
    """
    # Given parameters
    n = 1000
    time_float_op = 5  # in nanoseconds
    time_int_op = 1    # in nanoseconds
    time_func_call = 15 # in nanoseconds

    # --- Calculation for the FFT-based Algorithm ---
    print("Step 1: Calculate the execution time for the original FFT-based algorithm.")
    
    # Cost = (Number of function calls * time per call) + (Number of float ops * time per float op)
    num_func_calls_fft = 2 * n
    num_float_ops_fft = 4 * n
    
    time_fft = num_func_calls_fft * time_func_call + num_float_ops_fft * time_float_op
    
    print("The total time is the sum of time for function calls and floating-point operations.")
    print(f"Equation: Time_FFT = (2 * n) * time_func_call + (4 * n) * time_float_op")
    print(f"Substituting values: Time_FFT = (2 * {n}) * {time_func_call} + (4 * {n}) * {time_float_op}")
    print(f"Calculation: Time_FFT = {num_func_calls_fft} * {time_func_call} + {num_float_ops_fft} * {time_float_op}")
    print(f"Intermediate result: Time_FFT = {num_func_calls_fft * time_func_call} + {num_float_ops_fft * time_float_op}")
    print(f"Final result: Time_FFT = {time_fft} ns\n")


    # --- Calculation for the Direct Integer-based Algorithm ---
    print("Step 2: Calculate the execution time for the direct integer-based algorithm.")
    
    # Cost = (Number of float ops for conversion * time per float op) + (Number of integer ops * time per int op)
    num_float_ops_direct = 2 * n
    num_int_ops_direct = 2 * n**2
    
    time_direct = num_float_ops_direct * time_float_op + num_int_ops_direct * time_int_op
    
    print("The total time is the sum of time for data conversion and integer convolution.")
    print(f"Equation: Time_Direct = (2 * n) * time_float_op + (2 * n^2) * time_int_op")
    print(f"Substituting values: Time_Direct = (2 * {n}) * {time_float_op} + (2 * {n}**2) * {time_int_op}")
    print(f"Calculation: Time_Direct = {num_float_ops_direct} * {time_float_op} + {num_int_ops_direct} * {time_int_op}")
    print(f"Intermediate result: Time_Direct = {num_float_ops_direct * time_float_op} + {num_int_ops_direct * time_int_op}")
    print(f"Final result: Time_Direct = {time_direct} ns\n")


    # --- Comparison and Conclusion ---
    print("Step 3: Compare the two results.")
    print(f"FFT-based algorithm time: {time_fft} ns")
    print(f"Direct integer-based algorithm time: {time_direct} ns")

    if time_fft < time_direct:
        print("\nConclusion: The original FFT-based algorithm is faster.")
        answer = "Y"
    else:
        print("\nConclusion: The direct integer-based algorithm is faster.")
        answer = "N"
        
    print(f"\nIs the original FFT-based algorithm already the fastest? {answer}")

solve()