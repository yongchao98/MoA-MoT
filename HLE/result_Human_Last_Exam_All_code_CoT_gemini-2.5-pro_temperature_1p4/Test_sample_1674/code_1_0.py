import numpy as np
from scipy.linalg import expm

def run_simulation():
    """
    Simulates light propagation through an optical system and its time-reversal.
    Demonstrates the effect of an ideal vs. a dichroic birefringent medium.
    """
    # 1. Define the initial state (horizontally polarized light)
    E_in = np.array([[1], 
                     [0]], dtype=complex)
    E_in_conjugate = np.conj(E_in)
    
    # 2. Create a random medium 'T' (e.g., a depolarizer)
    # We ensure it's unitary (lossless) by creating a random Hermitian matrix H
    # and then calculating T = exp(iH).
    Z = np.random.randn(2, 2) + 1j * np.random.randn(2, 2)
    H = (Z + Z.conj().T) / 2  # Make it Hermitian
    T = expm(1j * H) # The unitary transmission matrix for the random medium

    # 3. Define the birefringent media
    # A quarter-wave plate (retardance of pi/2) is used as an example
    retardance = np.pi / 2
    
    # Case A: Ideal Birefringent Medium (unitary, no loss)
    B_ideal = np.array([[np.exp(-1j * retardance / 2), 0], 
                        [0, np.exp(1j * retardance / 2)]])

    # Case B: Birefringent Medium with Dichroism (non-unitary, has loss)
    # It has 100% transmission on x-axis and 50% on y-axis (p=0.5)
    p = 0.5 
    B_dichroic = np.array([[1 * np.exp(-1j * retardance / 2), 0], 
                           [0, p * np.exp(1j * retardance / 2)]])

    # --- Simulation 1: System with IDEAL Birefringent Medium ---
    print("--- Simulating Ideal System (Birefringence only) ---")
    # Total forward system matrix
    M_ideal = B_ideal @ T
    # Propagate forward
    E_out_ideal = M_ideal @ E_in
    # The "physical" time-reversal process
    E_out_ideal_conjugate = np.conj(E_out_ideal)
    M_ideal_reversed = M_ideal.T  # For a reciprocal medium, reverse operator is the transpose
    E_final_ideal = M_ideal_reversed @ E_out_ideal_conjugate

    print(f"Original Input (Conjugate): {E_in_conjugate.flatten()}")
    print(f"Reconstructed Input:        {E_final_ideal.flatten()}")
    # Check if the theory holds (vectors are numerically close)
    theory_holds_ideal = np.allclose(E_final_ideal, E_in_conjugate)
    print(f"Theory Holds: {theory_holds_ideal}\n")


    # --- Simulation 2: System with DICHROIC Birefringent Medium ---
    print("--- Simulating Non-Ideal System (Birefringence + Dichroism) ---")
    # Total forward system matrix
    M_dichroic = B_dichroic @ T
    # Propagate forward
    E_out_dichroic = M_dichroic @ E_in
    # The "physical" time-reversal process
    E_out_dichroic_conjugate = np.conj(E_out_dichroic)
    M_dichroic_reversed = M_dichroic.T
    E_final_dichroic = M_dichroic_reversed @ E_out_dichroic_conjugate

    print(f"Original Input (Conjugate): {E_in_conjugate.flatten()}")
    print(f"Reconstructed Input:        {E_final_dichroic.flatten()}")
    # Check if the theory holds
    theory_holds_dichroic = np.allclose(E_final_dichroic, E_in_conjugate)
    print(f"Theory Holds: {theory_holds_dichroic}")

if __name__ == "__main__":
    run_simulation()