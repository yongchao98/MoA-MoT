import numpy as np

def solve_dft():
    """
    Calculates the 8-point DFT of an interleaved sequence using the
    decimation-in-time FFT formulas.
    """
    # Given 4-point DFTs
    X = np.array([1, 1j, -1, -1j])
    H = np.array([0, 1 + 1j, 1, 1 - 1j])

    N = 8
    half_N = N // 2
    
    # Initialize the output 8-point DFT array
    Y = np.zeros(N, dtype=complex)

    print("This problem uses the decimation-in-time FFT structure.")
    print("The 8-point DFT Y(k) is calculated from the 4-point DFTs X(k) and H(k).")
    print("Formulas:")
    print("Y(k)   = X(k) + W_8^k * H(k)")
    print("Y(k+4) = X(k) - W_8^k * H(k)\n")
    print("--- Detailed Calculations ---\n")

    # Iterate from k = 0 to 3
    for k in range(half_N):
        # Calculate the twiddle factor W_8^k
        twiddle_factor = np.exp(-2j * np.pi * k / N)
        
        # Calculate the product term for clarity
        product_term = twiddle_factor * H[k]
        
        # --- Calculate Y(k) ---
        Y[k] = X[k] + product_term
        
        print(f"For k = {k}:")
        print(f"Y({k}) = X({k}) + W_8^{k} * H({k})")
        print(f"   = {X[k]} + ({twiddle_factor.real:.4f} + {twiddle_factor.imag:.4f}j) * {H[k]}")
        print(f"   = {X[k]} + ({product_term.real:.4f} + {product_term.imag:.4f}j)")
        print(f"   = {Y[k].real:.4f} + {Y[k].imag:.4f}j\n")

        # --- Calculate Y(k+4) ---
        Y[k + half_N] = X[k] - product_term
        
        print(f"For k = {k} (calculating Y({k+half_N})):")
        print(f"Y({k + half_N}) = X({k}) - W_8^{k} * H({k})")
        print(f"   = {X[k]} - ({twiddle_factor.real:.4f} + {twiddle_factor.imag:.4f}j) * {H[k]}")
        print(f"   = {X[k]} - ({product_term.real:.4f} + {product_term.imag:.4f}j)")
        print(f"   = {Y[k + half_N].real:.4f} + {Y[k + half_N].imag:.4f}j\n")
        print("--------------------------------\n")

    print("Final 8-point DFT Y(k) vector:")
    # Clean up near-zero floating point artifacts for printing
    Y_cleaned = [np.round(val, 10) for val in Y]
    print(np.array(Y_cleaned))
    
    # The final answer block as requested
    global final_answer
    final_answer = np.array(Y_cleaned)


# Run the solver
solve_dft()
# The final answer is wrapped in the specified format
# We convert complex numbers to strings for the final output format.
answer_str = "[" + ", ".join([f"({c.real}{c.imag:+}j)" for c in final_answer]) + "]"
print(f"<<<{answer_str}>>>")