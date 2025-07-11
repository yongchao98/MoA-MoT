import math

def calculate_convolution_times():
    """
    Calculates and compares the estimated execution time for different
    convolution algorithms on a specific machine.
    """
    # Machine operation times in nanoseconds (ns)
    time_int_add = 1
    time_int_mul = 2
    time_float_add = 9
    time_float_mul = 19

    # Problem size
    n = 1000

    print("Analysis of Convolution Algorithms for n = 1000\n")

    # --- Direct Convolution ---
    print("--- Method: Direct Convolution ---")
    num_muls_direct = n * n
    num_adds_direct = n * (n - 1)
    print(f"Operation Count: {num_muls_direct} multiplications, {num_adds_direct} additions")

    # Case B: Direct convolution with integers
    total_time_direct_int = (num_muls_direct * time_int_mul) + (num_adds_direct * time_int_add)
    print("\nB. Direct Convolution (Integers):")
    print(f"Equation: {num_muls_direct} * {time_int_mul} ns + {num_adds_direct} * {time_int_add} ns")
    print(f"Total Estimated Time = {total_time_direct_int:,} ns")

    # Case C: Direct convolution with floating points
    total_time_direct_float = (num_muls_direct * time_float_mul) + (num_adds_direct * time_float_add)
    print("\nC. Direct Convolution (Floating Points):")
    print(f"Equation: {num_muls_direct} * {time_float_mul} ns + {num_adds_direct} * {time_float_add} ns")
    print(f"Total Estimated Time = {total_time_direct_float:,} ns")

    # --- FFT-based Convolution ---
    print("\n--- Method: FFT-based Convolution ---")
    # FFT size N must be a power of 2 and >= 2n - 1
    output_len = 2 * n - 1
    N = 2**(math.ceil(math.log2(output_len)))
    log2_N = int(math.log2(N))
    print(f"FFT size N = {N}, log2(N) = {log2_N}")

    # Total complex operations for 3 FFTs and N element-wise products
    num_complex_muls = 3 * (N / 2) * log2_N + N
    num_complex_adds = 3 * N * log2_N

    # Convert complex operations to real floating-point operations
    # 1 complex mul = 4 real muls + 2 real adds
    # 1 complex add = 2 real adds
    num_float_muls_fft = int(num_complex_muls * 4)
    num_float_adds_fft = int((num_complex_muls * 2) + (num_complex_adds * 2))
    
    print(f"Operation Count: {num_float_muls_fft} float multiplications, {num_float_adds_fft} float additions")

    # Case A: FFT-based convolution (uses floating points)
    total_time_fft = (num_float_muls_fft * time_float_mul) + (num_float_adds_fft * time_float_add)
    print("\nA. FFT-based Convolution (Floating Points):")
    print(f"Equation: {num_float_muls_fft} * {time_float_mul} ns + {num_float_adds_fft} * {time_float_add} ns")
    print(f"Total Estimated Time = {total_time_fft:,} ns")

    # --- Conclusion ---
    print("\n--- Conclusion ---")
    times = {
        "Direct Convolution (Integers)": total_time_direct_int,
        "FFT-based Convolution": total_time_fft,
        "Direct Convolution (Floating Points)": total_time_direct_float,
    }
    fastest_method = min(times, key=times.get)
    print(f"The fastest algorithm is: {fastest_method}")
    print(f"Estimated time: {times[fastest_method]:,} ns")

if __name__ == '__main__':
    calculate_convolution_times()