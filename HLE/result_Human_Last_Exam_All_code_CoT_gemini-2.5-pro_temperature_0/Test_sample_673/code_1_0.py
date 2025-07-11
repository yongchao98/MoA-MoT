import math

def solve():
    """
    Compares the execution time of an FFT-based convolution algorithm
    with a direct integer-based convolution algorithm.
    """
    # --- Problem Parameters ---
    n = 1000  # Typical vector size
    T_float_op = 5  # ns
    T_int_op = 1    # ns
    T_call = 15     # ns

    # --- Algorithm 1: FFT-based ---
    # The time complexity is derived from the recursive structure. A recursive FFT of size n
    # makes 2 calls on size n/2. This happens log2(n) times down the recursion tree.
    # The problem states 4n floating point operations are performed at a step, which we interpret
    # as the work done at each level of the recursion.
    # Total float ops = (ops per element at one level) * n * (number of levels)
    # Total float ops = 4 * n * log2(n)
    # The number of function calls in a binary recursion on n items is 2*(n-1).
    num_float_ops_fft = 4 * n * math.log2(n)
    num_calls_fft = 2 * (n - 1)

    # Total time for FFT-based algorithm
    time_fft = num_float_ops_fft * T_float_op + num_calls_fft * T_call

    # --- Algorithm 2: Integer-based Direct Convolution ---
    # The problem states the costs directly: 2n float ops and 2n^2 int ops.
    num_float_ops_int = 2 * n
    num_int_ops_int = 2 * n**2

    # Total time for integer-based algorithm
    time_int = num_float_ops_int * T_float_op + num_int_ops_int * T_int_op

    # --- Output Results ---
    print(f"Comparing two algorithms for convolution of vectors of size n = {n}.")
    print("-" * 70)

    print("Algorithm 1: FFT-based method")
    print(f"Time(n) = (4 * n * log2(n)) * T_float_op + (2 * (n - 1)) * T_call")
    print(f"Time({n}) = (4 * {n} * log2({n})) * {T_float_op} ns + (2 * ({n} - 1)) * {T_call} ns")
    print(f"Time({n}) = (4 * {n} * {math.log2(n):.4f}) * {T_float_op} + (2 * {n-1}) * {T_call}")
    print(f"Time({n}) = {num_float_ops_fft * T_float_op:.2f} ns (float ops) + {num_calls_fft * T_call:.2f} ns (calls)")
    print(f"Total time for FFT-based method: {time_fft:.2f} ns")
    print("-" * 70)

    print("Algorithm 2: Integer-based direct convolution")
    print(f"Time(n) = (2 * n) * T_float_op + (2 * n^2) * T_int_op")
    print(f"Time({n}) = (2 * {n}) * {T_float_op} ns + (2 * {n}^2) * {T_int_op} ns")
    print(f"Time({n}) = ({num_float_ops_int}) * {T_float_op} + ({num_int_ops_int}) * {T_int_op}")
    print(f"Time({n}) = {num_float_ops_int * T_float_op:.2f} ns (float ops) + {num_int_ops_int * T_int_op:.2f} ns (int ops)")
    print(f"Total time for Integer-based method: {time_int:.2f} ns")
    print("-" * 70)

    print("Comparison:")
    if time_fft < time_int:
        print(f"The FFT-based algorithm ({time_fft:.2f} ns) is faster than the integer-based one ({time_int:.2f} ns).")
        answer = "Y"
    else:
        print(f"The integer-based algorithm ({time_int:.2f} ns) is faster than the FFT-based one ({time_fft:.2f} ns).")
        answer = "N"

    print("\nIs the original FFT-based algorithm already the fastest?")
    print(f"<<<{answer}>>>")

solve()