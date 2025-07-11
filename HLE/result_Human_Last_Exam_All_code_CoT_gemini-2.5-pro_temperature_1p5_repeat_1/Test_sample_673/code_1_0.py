import math

def solve():
    """
    Analyzes and compares the runtime of two convolution algorithms.
    """
    # Define the parameters from the problem description
    n = 1000
    float_op_time = 5  # ns
    int_op_time = 1    # ns
    func_call_time = 15 # ns

    # --- Algorithm 1: FFT-based Method ---
    # The algorithm is a divide-and-conquer type, characteristic of FFT.
    # The recurrence relation is T(n) = 2*T(n/2) + Work(n).
    # Work at each step of size 'm' consists of 4m float operations and 1 function call.
    # Time for this work is: 4 * m * float_op_time + func_call_time = 20m + 15.
    # The total time is the sum of work over all log2(n) recursion levels.
    # At recursion level k (from 0), there are 2^k subproblems, each of size n/2^k.
    # Work at level k = 2^k * (20 * (n/2^k) + 15) = 20*n + 15*2^k.
    # Total time = Sum_{k=0 to log2(n)-1} (20*n + 15*2^k)
    # This simplifies to: 20*n*log2(n) + 15*(n-1).

    log2_n = math.log2(n)
    cost_fft = (20 * n * log2_n) + (15 * (n - 1))

    # --- Algorithm 2: Direct Integer Convolution ---
    # This method involves a sequence of operations:
    # 1. Conversion: 2n floating point operations.
    # 2. Convolution: 2n^2 integer operations.
    # 3. Overhead: We assume a single function call for the entire process.
    # Total time = (2*n * float_op_time) + (2*n^2 * int_op_time) + func_call_time

    cost_direct = (2 * n * float_op_time) + (2 * n**2 * int_op_time) + func_call_time
    
    # --- Output the results ---
    print("Runtime Analysis of Convolution Algorithms")
    print("==========================================")
    print(f"Given Parameters:")
    print(f"  Vector size, n = {n}")
    print(f"  Float operation time = {float_op_time} ns")
    print(f"  Integer operation time = {int_op_time} ns")
    print(f"  Function call time = {func_call_time} ns\n")

    print("--- Algorithm 1: FFT-based Algorithm ---")
    print("Formula: Time = (20 * n * log2(n)) + (15 * (n - 1))")
    print(f"Calculation: Time = (20 * {n} * {log2_n:.4f}) + (15 * ({n} - 1))")
    print(f"Total Estimated Time: {cost_fft:,.2f} ns\n")
    
    print("--- Algorithm 2: Direct Integer Convolution ---")
    print("Formula: Time = (2 * n * {float_op_time}) + (2 * n^2 * {int_op_time}) + {func_call_time}")
    print(f"Calculation: Time = (2 * {n} * {float_op_time}) + (2 * {n}*{n} * {int_op_time}) + {func_call_time}")
    print(f"Total Estimated Time: {cost_direct:,.2f} ns\n")

    print("--- Conclusion ---")
    if cost_fft < cost_direct:
        print(f"The FFT-based algorithm ({cost_fft:,.2f} ns) is faster than the direct algorithm ({cost_direct:,.2f} ns).")
        answer = "Y"
    else:
        print(f"The direct algorithm ({cost_direct:,.2f} ns) is faster than the FFT-based algorithm ({cost_fft:,.2f} ns).")
        answer = "N"
        
    print(f"\nIs the original FFT-based algorithm already the fastest?")
    print(f"<<<{answer}>>>")

solve()