import math

def solve():
    # --- Problem Definition ---
    n = 1000
    t_iadd = 1  # ns
    t_imul = 2  # ns
    t_fadd = 9  # ns
    t_fmul = 19 # ns

    # --- Introduction ---
    print("Analysis of Convolution Algorithms")
    print("-" * 40)
    print("This script estimates the execution time for convoluting two series of n=1000 elements on a specific machine.")
    print("\nMachine Operation Times:")
    print(f"- Integer Addition/Subtraction: {t_iadd} ns")
    print(f"- Integer Multiplication: {t_imul} ns")
    print(f"- Floating Point Addition/Subtraction: {t_fadd} ns")
    print(f"- Floating Point Multiplication: {t_fmul} ns")
    print(f"\nSeries size n = {n} elements.")
    print("-" * 40)

    # --- Method 1: Direct Convolution ---
    print("\n--- METHOD 1: Direct Convolution ---")
    print("This method computes the convolution sum y[m] = sum(x[k] * h[m-k]) directly.")
    print(f"For two series of length n={n}, this requires approximately n^2 multiplications and n^2 additions.")

    n_ops_direct_mul = n ** 2
    n_ops_direct_add = n ** 2

    print("\nNumber of operations:")
    print(f"- Multiplications (n^2): {n}^2 = {n_ops_direct_mul}")
    print(f"- Additions (n^2): {n}^2 = {n_ops_direct_add}")

    # B. Direct convolution with integers
    time_direct_int_mul = n_ops_direct_mul * t_imul
    time_direct_int_add = n_ops_direct_add * t_iadd
    total_time_direct_int = time_direct_int_mul + time_direct_int_add

    print("\nB. Direct convolution with integers:")
    print(f"- Time for multiplications = {n_ops_direct_mul} * {t_imul} ns = {time_direct_int_mul} ns")
    print(f"- Time for additions       = {n_ops_direct_add} * {t_iadd} ns = {time_direct_int_add} ns")
    print(f"- Total time = {time_direct_int_mul} + {time_direct_int_add} = {total_time_direct_int} ns ({total_time_direct_int / 1e6:.2f} ms)")

    # C. Direct convolution with floating points
    time_direct_float_mul = n_ops_direct_mul * t_fmul
    time_direct_float_add = n_ops_direct_add * t_fadd
    total_time_direct_float = time_direct_float_mul + time_direct_float_add

    print("\nC. Direct convolution with floating points:")
    print(f"- Time for multiplications = {n_ops_direct_mul} * {t_fmul} ns = {time_direct_float_mul} ns")
    print(f"- Time for additions       = {n_ops_direct_add} * {t_fadd} ns = {time_direct_float_add} ns")
    print(f"- Total time = {time_direct_float_mul} + {time_direct_float_add} = {total_time_direct_float} ns ({total_time_direct_float / 1e6:.2f} ms)")
    print("-" * 40)

    # --- Method 2: FFT-based Convolution ---
    print("\n--- METHOD 2: FFT-based Convolution (float) ---")
    print("This method uses the Convolution Theorem. Steps: FFT(x), FFT(h), Pointwise-Multiply, IFFT(result).")
    
    N_padded = 2 * n - 1
    N_fft = 1
    while N_fft < N_padded:
        N_fft *= 2
    log2N = int(math.log2(N_fft))

    print("\nFFT Parameters:")
    print(f"- Required length for full convolution: 2*n - 1 = {N_padded}")
    print(f"- Next power of 2 for FFT efficiency, N = {N_fft}")
    print(f"- log2(N) = {log2N}")
    
    print("\nWe assume an optimized algorithm for real-valued signals.")
    # For a complex N-point FFT, ops are approx:
    # Muls: 2*N*log2(N), Adds: 3*N*log2(N)
    # A real N-point FFT takes about half the operations.
    muls_per_real_fft = (1/2) * (2 * N_fft * log2N) # approx
    adds_per_real_fft = (1/2) * (3 * N_fft * log2N) # approx

    # Pointwise multiplication of two complex vectors with Hermitian symmetry
    # Y[k] = X[k] * H[k] for k=0..N/2
    muls_pointwise = 4 * (N_fft // 2 - 1) + 2
    adds_pointwise = 2 * (N_fft // 2 - 1)

    # Total ops = 2 forward real FFTs + 1 inverse real FFT + 1 pointwise multiply
    total_fft_muls = 2 * muls_per_real_fft + muls_pointwise + muls_per_real_fft
    total_fft_adds = 2 * adds_per_real_fft + adds_pointwise + adds_per_real_fft
    
    # Casting to int for cleaner output
    total_fft_muls = int(total_fft_muls)
    total_fft_adds = int(total_fft_adds)
    
    # A. FFT-based convolution with floating points
    time_fft_mul = total_fft_muls * t_fmul
    time_fft_add = total_fft_adds * t_fadd
    total_time_fft = time_fft_mul + time_fft_add

    print("\nA. FFT-based convolution with floating points:")
    print("Total operations comprise 2 forward FFTs, 1 inverse FFT, and 1 point-wise multiplication.")
    print(f"- Total Float Multiplications = {total_fft_muls}")
    print(f"- Total Float Additions       = {total_fft_adds}")
    print(f"- Time for multiplications = {total_fft_muls} * {t_fmul} ns = {time_fft_mul} ns")
    print(f"- Time for additions       = {total_fft_adds} * {t_fadd} ns = {time_fft_add} ns")
    print(f"- Total time = {time_fft_mul} + {time_fft_add} = {total_time_fft} ns ({total_time_fft / 1e6:.2f} ms)")
    print("-" * 40)

    # --- Conclusion ---
    print("\n--- Conclusion ---")
    print("Comparing the total estimated times:")
    print(f"- A. FFT (float):               {total_time_fft / 1e6:.2f} ms")
    print(f"- B. Direct (integer):          {total_time_direct_int / 1e6:.2f} ms")
    print(f"- C. Direct (float):            {total_time_direct_float / 1e6:.2f} ms")

    print("\nThe fastest algorithm is FFT-based convolution. Its O(N log N) complexity is more")
    print("efficient than the O(n^2) direct method for a large series like n=1000, even though")
    print("it uses more expensive floating-point operations.")

if __name__ == '__main__':
    solve()
<<<A>>>