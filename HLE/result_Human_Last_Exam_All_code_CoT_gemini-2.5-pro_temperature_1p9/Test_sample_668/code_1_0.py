import math

def solve_convolution_time():
    """
    Calculates and compares the execution time for different convolution algorithms.
    """
    # Machine operation times in nanoseconds (ns)
    t_iadd = 1
    t_imul = 2
    t_fadd = 9
    t_fmul = 19
    
    # Parameters of the input series
    n = 1000

    print("Step-by-step calculation to find the fastest convolution algorithm.")
    print("The machine has the following operation times (in nanoseconds):")
    print(f"- integer addition/subtraction: {t_iadd} ns")
    print(f"- integer multiplication: {t_imul} ns")
    print(f"- floating point addition/subtraction: {t_fadd} ns")
    print(f"- floating point multiplication: {t_fmul} ns")
    print(f"We need to find the fastest method to convolve two series of n = {n} elements.")
    print("")

    # --- Analysis 1: Direct Convolution ---
    print("--- Method 1: Direct Convolution (O(n^2) complexity) ---")
    print(f"Direct convolution of two n-element series requires approximately n^2 multiplications and n^2 additions.")
    num_ops_direct = n**2
    print(f"For n = {n}, we need {n}*{n} = {num_ops_direct} operations of each type.")
    print("")

    print("A. Calculation for Direct Convolution with Integers (Choice B):")
    time_mul_int = num_ops_direct * t_imul
    time_add_int = num_ops_direct * t_iadd
    time_direct_int = time_mul_int + time_add_int
    print(f"Total Time = ( {num_ops_direct} multiplications * {t_imul} ns/mul ) + ( {num_ops_direct} additions * {t_iadd} ns/add )")
    print(f"Total Time = {time_mul_int} ns + {time_add_int} ns = {time_direct_int} ns")
    print("")

    print("B. Calculation for Direct Convolution with Floating Points (Choice C):")
    time_mul_float = num_ops_direct * t_fmul
    time_add_float = num_ops_direct * t_fadd
    time_direct_float = time_mul_float + time_add_float
    print(f"Total Time = ( {num_ops_direct} multiplications * {t_fmul} ns/mul ) + ( {num_ops_direct} additions * {t_fadd} ns/add )")
    print(f"Total Time = {time_mul_float} ns + {time_add_float} ns = {time_direct_float} ns")
    print("")

    # --- Analysis 2: FFT-based Convolution ---
    print("--- Method 2: FFT-based Convolution (O(N log N) complexity) (Choice A) ---")
    print("This method uses the Convolution Theorem: convolution(x, h) = IFFT( FFT(x) * FFT(h) ).")
    print("It involves floating point arithmetic with complex numbers.")

    # Step 2.1: Find the padded length N
    output_len = 2 * n - 1
    N = 1
    while N < output_len:
        N *= 2
    log2_N = int(math.log2(N))

    print(f"1. Pad sequences: The sequences must be padded to a length N >= 2*n-1 = {output_len}. N must be a power of 2.")
    print(f"The chosen padded length is N = {N}.")
    print("")

    print("2. Operation Count: An N-point FFT/IFFT uses about (N/2)*log2(N) complex multiplications and N*log2(N) complex additions.")
    # A single complex multiplication (a+ib)*(c+id) needs 4 real muls and 2 real adds.
    # A single complex addition needs 2 real adds.
    num_cmul_fft = (N / 2) * log2_N
    num_cadd_fft = N * log2_N
    muls_per_fft = num_cmul_fft * 4
    adds_per_fft = (num_cmul_fft * 2) + (num_cadd_fft * 2)
    muls_per_ew_prod = N * 4
    adds_per_ew_prod = N * 2

    print("Cost of sub-steps in terms of real floating-point operations:")
    print(f" - For one {N}-point FFT or IFFT:")
    print(f"   - Complex Multiplications: ({N}/2)*{log2_N} = {int(num_cmul_fft)}")
    print(f"   - Complex Additions: {N}*{log2_N} = {int(num_cadd_fft)}")
    print(f"   - Equivalent Real Multiplications: {int(num_cmul_fft)} * 4 = {int(muls_per_fft)}")
    print(f"   - Equivalent Real Additions: ({int(num_cmul_fft)} * 2) + ({int(num_cadd_fft)} * 2) = {int(adds_per_fft)}")
    print(f" - For the element-wise product of {N} complex numbers:")
    print(f"   - Equivalent Real Multiplications: {N} * 4 = {int(muls_per_ew_prod)}")
    print(f"   - Equivalent Real Additions: {N} * 2 = {int(adds_per_ew_prod)}")
    print("")

    print("3. Total Operation Count for the full FFT convolution (FFT(x) + FFT(h) + Product + IFFT(Y)):")
    # The process requires 2 FFTs, 1 IFFT, and 1 element-wise product.
    total_muls_fft = muls_per_fft + muls_per_fft + muls_per_ew_prod + muls_per_fft
    total_adds_fft = adds_per_fft + adds_per_fft + adds_per_ew_prod + adds_per_fft
    print(f"Total Real Multiplications = Muls(FFT_x) + Muls(FFT_h) + Muls(Product) + Muls(IFFT_y)")
    print(f"Total Real Multiplications = {int(muls_per_fft)} + {int(muls_per_fft)} + {int(muls_per_ew_prod)} + {int(muls_per_fft)} = {int(total_muls_fft)}")
    print(f"Total Real Additions = Adds(FFT_x) + Adds(FFT_h) + Adds(Product) + Adds(IFFT_y)")
    print(f"Total Real Additions = {int(adds_per_fft)} + {int(adds_per_fft)} + {int(adds_per_ew_prod)} + {int(adds_per_fft)} = {int(total_adds_fft)}")
    print("")

    print("4. Final time calculation for FFT-based convolution:")
    time_fft_mul = total_muls_fft * t_fmul
    time_fft_add = total_adds_fft * t_fadd
    time_fft_float = time_fft_mul + time_fft_add
    print(f"Total Time = ( {int(total_muls_fft)} multiplications * {t_fmul} ns/mul ) + ( {int(total_adds_fft)} additions * {t_fadd} ns/add )")
    print(f"Total Time = {int(time_fft_mul)} ns + {int(time_fft_add)} ns = {int(time_fft_float)} ns")
    print("")

    # --- Conclusion ---
    print("--- Conclusion ---")
    print("Comparing the total estimated times:")
    print(f"- Direct Convolution (Integers):     {time_direct_int:,} ns")
    print(f"- FFT Convolution (Floating Points): {int(time_fft_float):,} ns")
    print(f"- Direct Convolution (Floating Points): {time_direct_float:,} ns")
    print("")

    if time_direct_int < time_fft_float and time_direct_int < time_direct_float:
        print("The fastest method is Direct convolution with integers.")
    elif time_fft_float < time_direct_int and time_fft_float < time_direct_float:
        print("The fastest method is FFT-based convolution.")
    else:
        print("The fastest method is Direct convolution with floating points.")

if __name__ == '__main__':
    solve_convolution_time()
