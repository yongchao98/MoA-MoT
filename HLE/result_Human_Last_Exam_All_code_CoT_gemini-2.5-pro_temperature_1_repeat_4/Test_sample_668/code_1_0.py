import math

def calculate_convolution_times():
    """
    Calculates and prints the estimated execution times for different
    convolution algorithms based on given operation costs.
    """

    # Machine operation times in nanoseconds
    T_INT_ADD = 1
    T_INT_MUL = 2
    T_FLOAT_ADD = 9
    T_FLOAT_MUL = 19

    # Series length
    n = 1000

    print("Step-by-step calculation of convolution time:")
    print("--------------------------------------------------")

    # --- Case 1: Direct Convolution with Integers ---
    print("1. Direct Convolution with Integers (Time in ns)")
    # Complexity: O(n^2) multiplications and O(n^2) additions.
    num_mult_direct = n**2
    num_add_direct = n**2
    time_direct_int = num_mult_direct * T_INT_MUL + num_add_direct * T_INT_ADD
    
    print(f"Number of integer multiplications = {n} * {n} = {num_mult_direct}")
    print(f"Number of integer additions = {n} * {n} = {num_add_direct}")
    print(f"Total time = ({num_mult_direct} * {T_INT_MUL}) + ({num_add_direct} * {T_INT_ADD}) = {num_mult_direct * T_INT_MUL + num_add_direct * T_INT_ADD}")
    print("--------------------------------------------------")

    # --- Case 2: Direct Convolution with Floating Points ---
    print("2. Direct Convolution with Floating Points (Time in ns)")
    # Same number of operations, but with floating point costs.
    time_direct_float = num_mult_direct * T_FLOAT_MUL + num_add_direct * T_FLOAT_ADD

    print(f"Number of floating point multiplications = {n} * {n} = {num_mult_direct}")
    print(f"Number of floating point additions = {n} * {n} = {num_add_direct}")
    print(f"Total time = ({num_mult_direct} * {T_FLOAT_MUL}) + ({num_add_direct} * {T_FLOAT_ADD}) = {num_mult_direct * T_FLOAT_MUL + num_add_direct * T_FLOAT_ADD}")
    print("--------------------------------------------------")

    # --- Case 3: FFT-based Convolution (with Floating Points) ---
    print("3. FFT-based Convolution (Time in ns)")
    # Padded length N must be a power of 2 >= 2n - 1.
    required_len = 2 * n - 1
    N = 1
    while N < required_len:
        N *= 2
    log2N = int(math.log2(N))

    # Ops for one FFT: (N/2)*log2(N) complex mults, N*log2(N) complex adds.
    # Total ops: 3 FFTs + N element-wise complex mults.
    total_complex_mults = 3 * (N / 2 * log2N) + N
    total_complex_adds = 3 * (N * log2N)

    # Convert complex ops to real ops:
    # 1 complex mult = 4 real mults + 2 real adds
    # 1 complex add = 2 real adds
    total_float_mults = total_complex_mults * 4
    total_float_adds = (total_complex_mults * 2) + (total_complex_adds * 2)
    time_fft = total_float_mults * T_FLOAT_MUL + total_float_adds * T_FLOAT_ADD

    print(f"Required transform length >= 2*{n} - 1 = {required_len}")
    print(f"Chosen FFT length N (next power of 2) = {N}")
    print(f"For this N, log2(N) = {log2N}")
    print("Total complex operations for 2 FFTs, 1 IFFT, and element-wise product:")
    print(f"Total Complex Multiplications = (3 * {N}/2 * {log2N}) + {N} = {int(total_complex_mults)}")
    print(f"Total Complex Additions = 3 * {N} * {log2N} = {int(total_complex_adds)}")
    print("Converting to real floating point operations:")
    print(f"Total Real Multiplications = {int(total_complex_mults)} * 4 = {int(total_float_mults)}")
    print(f"Total Real Additions = ({int(total_complex_mults)} * 2) + ({int(total_complex_adds)} * 2) = {int(total_float_adds)}")
    print("Final Time Calculation:")
    print(f"Total time = ({int(total_float_mults)} * {T_FLOAT_MUL}) + ({int(total_float_adds)} * {T_FLOAT_ADD}) = {int(time_fft)}")
    print("--------------------------------------------------")

    # --- Comparison ---
    print("Comparison of Total Times:")
    print(f"Direct Integer Convolution:          {time_direct_int} ns")
    print(f"Direct Floating-Point Convolution:   {time_direct_float} ns")
    print(f"FFT-based Convolution:               {int(time_fft)} ns")
    print("--------------------------------------------------")

    if time_direct_int < time_direct_float and time_direct_int < time_fft:
        fastest = "Direct convolution with integers"
    elif time_direct_float < time_fft:
        fastest = "Direct convolution with floating points"
    else:
        fastest = "FFT-based convolution"

    print(f"Conclusion: The fastest algorithm is '{fastest}'.")

# Execute the calculation
calculate_convolution_times()