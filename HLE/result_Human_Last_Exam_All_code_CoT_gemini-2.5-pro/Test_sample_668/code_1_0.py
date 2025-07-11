import math

def solve_convolution_time():
    # Machine operation times in nanoseconds (ns)
    time_int_add = 1
    time_int_mul = 2
    time_fp_add = 9
    time_fp_mul = 19

    # Series size
    n = 1000

    print("Analysis of Convolution Time for n = 1000")
    print("==========================================")
    print("Machine Operation Times:")
    print(f"- Integer Addition/Subtraction: {time_int_add} ns")
    print(f"- Integer Multiplication: {time_int_mul} ns")
    print(f"- Floating Point Addition/Subtraction: {time_fp_add} ns")
    print(f"- Floating Point Multiplication: {time_fp_mul} ns")
    print("\n--- Method 1: Direct Convolution ---\n")

    # --- Direct Convolution ---
    # Number of operations is approximately n^2 for both additions and multiplications.
    direct_ops_count = n * n

    # 1.A. Direct convolution with integers
    print("A. Direct Convolution with Integers")
    time_direct_int = (direct_ops_count * time_int_mul) + (direct_ops_count * time_int_add)
    print(f"Number of multiplications = {n} * {n} = {direct_ops_count}")
    print(f"Number of additions = {n} * {n} = {direct_ops_count}")
    print(f"Total Time = ({direct_ops_count} * {time_int_mul} ns) + ({direct_ops_count} * {time_int_add} ns) = {time_direct_int:,} ns")
    print("\n")

    # 1.B. Direct convolution with floating points
    print("B. Direct Convolution with Floating Points")
    time_direct_fp = (direct_ops_count * time_fp_mul) + (direct_ops_count * time_fp_add)
    print(f"Number of multiplications = {n} * {n} = {direct_ops_count}")
    print(f"Number of additions = {n} * {n} = {direct_ops_count}")
    print(f"Total Time = ({direct_ops_count} * {time_fp_mul} ns) + ({direct_ops_count} * {time_fp_add} ns) = {time_direct_fp:,} ns")
    print("\n")

    # --- FFT-based Convolution ---
    print("--- Method 2: FFT-based Convolution ---\n")
    # For linear convolution, pad to N >= 2n - 1
    min_N = 2 * n - 1
    # Find the next power of 2
    N = 2**math.ceil(math.log2(min_N))
    log2_N = int(math.log2(N))
    print(f"To avoid circular convolution, we pad the series from n={n} to N={N} (the next power of 2 >= {min_N}).")
    print(f"log2({N}) = {log2_N}")

    # Operations for 3 FFTs (2 forward, 1 inverse) + element-wise product
    # Number of complex multiplications for one FFT = (N/2)*log2(N)
    # Number of complex additions for one FFT = N*log2(N)
    
    # Total complex operations
    fft_complex_mul = 3 * (N / 2) * log2_N
    fft_complex_add = 3 * N * log2_N
    elementwise_complex_mul = N
    
    total_complex_mul = fft_complex_mul + elementwise_complex_mul
    total_complex_add = fft_complex_add
    
    print(f"\nTotal operations for 3 FFTs and 1 element-wise product:")
    print(f"Total complex multiplications = 3 * ({N}/2)*{log2_N} + {N} = {int(total_complex_mul)}")
    print(f"Total complex additions = 3 * {N}*{log2_N} = {int(total_complex_add)}")

    # Convert complex operations to real floating-point operations
    # 1 complex mul = 4 real mul + 2 real add
    # 1 complex add = 2 real add
    total_real_mul = total_complex_mul * 4
    total_real_add = (total_complex_mul * 2) + (total_complex_add * 2)
    
    print("\nConverting to real floating-point operations:")
    print(f"Total real multiplications = {int(total_complex_mul)} * 4 = {int(total_real_mul)}")
    print(f"Total real additions = ({int(total_complex_mul)} * 2) + ({int(total_complex_add)} * 2) = {int(total_real_add)}")

    # Calculate final time for FFT method
    time_fft = (total_real_mul * time_fp_mul) + (total_real_add * time_fp_add)
    print("\nC. Total Time for FFT-based Convolution")
    print(f"Total Time = ({int(total_real_mul)} * {time_fp_mul} ns) + ({int(total_real_add)} * {time_fp_add} ns) = {int(time_fft):,} ns")
    print("\n")

    # --- Conclusion ---
    print("--- Conclusion ---")
    print("\nSummary of Times:")
    print(f"A. FFT: {int(time_fft):,} ns")
    print(f"B. Direct convolution with integers: {time_direct_int:,} ns")
    print(f"C. Direct convolution with floating points: {time_direct_fp:,} ns")

    times = {
        "FFT": time_fft,
        "Direct convolution with integers": time_direct_int,
        "Direct convolution with floating points": time_direct_fp,
    }

    fastest_method = min(times, key=times.get)
    fastest_time = times[fastest_method]

    print(f"\nThe fastest algorithm is '{fastest_method}' with an estimated time of {fastest_time:,} ns.")

solve_convolution_time()
<<<B>>>