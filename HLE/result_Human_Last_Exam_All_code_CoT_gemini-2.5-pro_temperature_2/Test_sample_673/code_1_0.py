import math

def solve():
    """
    Compares the execution time of an FFT-based convolution algorithm
    with a direct integer-based convolution algorithm.
    """
    # Given parameters
    n = 1000
    t_float = 5  # ns
    t_int = 1    # ns
    t_call = 15  # ns

    # --- Algorithm 1: FFT-based ---
    # Time = (Function call overhead) + (Computational cost)
    # Function call overhead for a recursive function splitting in two is (2n-1) calls.
    # Computational cost is 4n*log2(n) floating point operations.
    
    t1_call_component = (2 * n - 1) * t_call
    # math.log2(1000) is approx 9.96578
    t1_ops_component = 4 * n * math.log2(n) * t_float
    t1_total = t1_call_component + t1_ops_component

    print("--- Algorithm 1: FFT-based ---")
    print(f"Time (T1) = (2 * n - 1) * t_call + 4 * n * log2(n) * t_float")
    print(f"T1 = (2 * {n} - 1) * {t_call} + 4 * {n} * {math.log2(n):.4f} * {t_float}")
    print(f"T1 = {t1_call_component:.0f} (call cost) + {t1_ops_component:.0f} (ops cost)")
    print(f"T1 = {t1_total:.0f} ns\n")

    # --- Algorithm 2: Integer-based direct convolution ---
    # Time = (Function call overhead) + (Conversion cost) + (Convolution cost)
    # A single function call.
    # Conversion: 2n float operations.
    # Convolution: 2n^2 integer operations.
    
    t2_call_component = 1 * t_call
    t2_conversion_component = 2 * n * t_float
    t2_convolution_component = 2 * n**2 * t_int
    t2_total = t2_call_component + t2_conversion_component + t2_convolution_component
    
    print("--- Algorithm 2: Integer-based ---")
    print(f"Time (T2) = t_call + 2 * n * t_float + 2 * n^2 * t_int")
    print(f"T2 = {t_call} + 2 * {n} * {t_float} + 2 * {n}^2 * {t_int}")
    print(f"T2 = {t2_call_component} (call cost) + {t2_conversion_component} (conversion cost) + {t2_convolution_component} (convolution cost)")
    print(f"T2 = {t2_total:.0f} ns\n")

    # --- Comparison ---
    is_original_faster = t1_total < t2_total
    
    print("--- Conclusion ---")
    if is_original_faster:
        print("The original FFT-based algorithm is faster.")
        answer = "Y"
    else:
        print("The integer-based direct convolution algorithm is faster.")
        answer = "N"
    
    # The final answer in the required format is printed outside this function.
    return answer

# Execute the analysis and store the final answer.
final_answer = solve()
print(f"\nIs the original FFT-based algorithm is already the fastest?\n{final_answer}")
