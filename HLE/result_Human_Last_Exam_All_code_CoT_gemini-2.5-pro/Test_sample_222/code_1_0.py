import numpy as np

def solve_dft():
    """
    Calculates the 4-point DFT of an 8-point sequence formed by interleaving
    two 4-point sequences, given their respective 4-point DFTs.
    """
    # Given 4-point DFTs
    X = np.array([1, 1j, -1, -1j])
    H = np.array([0, 1 + 1j, 1, 1 - 1j])

    print("This script calculates the 4-point DFT (let's call it Y4) of an 8-point sequence.")
    print("The 8-point sequence is formed by interleaving two 4-point sequences x(n) and h(n).")
    print("We are given the 4-point DFTs X(k) and H(k).")
    print("\nThe derived formula is: Y4(k) = X(2k mod 4) + W_4^k * H(2k mod 4) for k=0,1,2,3\n")

    # Array to store the results
    Y4 = np.zeros(4, dtype=complex)

    # Loop through k = 0, 1, 2, 3
    for k in range(4):
        # Twiddle factor W_4^k
        W4_k = np.exp(-1j * 2 * np.pi * k / 4)
        
        # Index for X and H is (2*k) mod 4 because they are 4-periodic
        idx = (2 * k) % 4
        
        X_val = X[idx]
        H_val = H[idx]
        
        # Calculate the result for the current k
        Y4[k] = X_val + W4_k * H_val
        
        # Format numbers for clean printing
        # Using round to avoid floating point artifacts like 1-0j
        x_str = f"({np.round(X_val, 4)})"
        h_str = f"({np.round(H_val, 4)})"
        w_str = f"({np.round(W4_k, 4)})"
        res_str = f"({np.round(Y4[k], 4)})"
        
        # Print the detailed calculation for the current step
        print(f"For k = {k}:")
        print(f"Y4({k}) = X({idx}) + W_4^{k} * H({idx})")
        print(f"      = {x_str} + {w_str} * {h_str}")
        print(f"      = {res_str}\n")

    print("The final result for the 4-point DFT is:")
    print(np.round(Y4, 4))

solve_dft()
<<<[ 1.+0.j -1.-1.j  1.+0.j -1.+1.j]>>>