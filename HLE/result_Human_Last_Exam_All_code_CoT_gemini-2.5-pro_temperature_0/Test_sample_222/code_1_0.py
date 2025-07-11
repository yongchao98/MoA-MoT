import numpy as np

def solve_dft_problem():
    """
    Calculates the 4-point DFT for the interleaved sequence as described.
    """
    # Given 4-point DFTs
    X = np.array([1, 1j, -1, -1j], dtype=complex)
    H = np.array([0, 1+1j, 1, 1-1j], dtype=complex)

    print("The problem is to find the 4-point DFT of the 8-point sequence y(n) formed by interleaving x(n) and h(n).")
    print("This is interpreted as finding the even-indexed components of the 8-point DFT of y(n), which we call Z(k).")
    print("Z(k) = {Y(0), Y(2), Y(4), Y(6)}\n")
    print("The calculation steps are as follows:\n")

    # Twiddle factors needed: W_8^0 = 1 and W_8^2 = -j
    W_8_0 = np.exp(-1j * 2 * np.pi * 0 / 8)
    W_8_2 = np.exp(-1j * 2 * np.pi * 2 / 8)

    # --- Calculate Z(0) = Y(0) ---
    k0 = 0
    # Y(0) = X(0) + W_8^0 * H(0)
    Y0 = X[k0] + W_8_0 * H[k0]
    print(f"Z(0) = Y(0) = X(0) + W_8^0 * H(0) = {X[k0]} + (1) * {H[k0]} = {np.round(Y0, 10)}")

    # --- Calculate Z(1) = Y(2) ---
    k2 = 2
    # Y(2) = X(2) + W_8^2 * H(2)
    Y2 = X[k2] + W_8_2 * H[k2]
    print(f"Z(1) = Y(2) = X(2) + W_8^2 * H(2) = {X[k2]} + (-1j) * {H[k2]} = {np.round(Y2, 10)}")

    # --- Calculate Z(2) = Y(4) ---
    # Y(4) = X(0) - W_8^0 * H(0)
    Y4 = X[k0] - W_8_0 * H[k0]
    print(f"Z(2) = Y(4) = X(0) - W_8^0 * H(0) = {X[k0]} - (1) * {H[k0]} = {np.round(Y4, 10)}")

    # --- Calculate Z(3) = Y(6) ---
    # Y(6) = X(2) - W_8^2 * H(2)
    Y6 = X[k2] - W_8_2 * H[k2]
    print(f"Z(3) = Y(6) = X(2) - W_8^2 * H(2) = {X[k2]} - (-1j) * {H[k2]} = {np.round(Y6, 10)}")

    # --- Final Result ---
    Z = [Y0, Y2, Y4, Y6]
    Z_cleaned = np.round(Z, 10)
    
    print("\n--------------------------------------------------")
    print("The final resulting 4-point DFT is:")
    print(f"Z = [{Z_cleaned[0]}, {Z_cleaned[1]}, {Z_cleaned[2]}, {Z_cleaned[3]}]")
    print("--------------------------------------------------")

solve_dft_problem()
<<<[1, (-1-1j), 1, (-1+1j)]>>>