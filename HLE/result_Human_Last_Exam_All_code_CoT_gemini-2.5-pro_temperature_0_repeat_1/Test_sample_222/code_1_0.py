import numpy as np

def solve_dft():
    """
    Calculates the 4-point DFT of an interleaved 8-point sequence.
    """
    # Given 4-point DFTs
    X = np.array([1, 1j, -1, -1j], dtype=complex)
    H = np.array([0, 1+1j, 1, 1-1j], dtype=complex)

    # The resulting 4-point DFT
    Z = np.zeros(4, dtype=complex)
    N = 4

    print("The 4-point DFT Z(k) of the interleaved sequence is calculated using the formula:")
    print("Z(k) = X((2*k) % 4) + exp(-j*2*pi*k/4) * H((2*k) % 4)\n")

    for k in range(N):
        # Twiddle factor W_4^k
        W = np.exp(-2j * np.pi * k / N)
        # Index for X and H, which is (2*k) mod 4
        idx = (2 * k) % N
        
        # Calculate Z(k)
        Z[k] = X[idx] + W * H[idx]
        
        # Format the numbers for clear output
        w_str = f"({W.real:.2f}{W.imag:+.2f}j)"
        x_val_str = f"({X[idx].real:.2f}{X[idx].imag:+.2f}j)"
        h_val_str = f"({H[idx].real:.2f}{H[idx].imag:+.2f}j)"
        
        # Clean up small floating point errors for the final result
        z_real = round(Z[k].real, 2)
        z_imag = round(Z[k].imag, 2)
        z_val_str = f"({z_real:.2f}{z_imag:+.2f}j)"
        
        print(f"For k = {k}:")
        print(f"  Z({k}) = X({idx}) + (W_4^{k}) * H({idx})")
        print(f"       = {x_val_str} + {w_str} * {h_val_str}")
        print(f"       = {z_val_str}\n")

    # Final result array
    final_Z = [complex(round(z.real, 2), round(z.imag, 2)) for z in Z]
    print("The final result for the 4-point DFT Z(k) is:")
    print(f"Z = {final_Z}")

solve_dft()
<<<[1+0j, -1-1j, 1+0j, -1+1j]>>>