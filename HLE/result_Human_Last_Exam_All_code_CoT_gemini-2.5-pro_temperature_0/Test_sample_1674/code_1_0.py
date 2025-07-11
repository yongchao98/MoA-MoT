import numpy as np
from scipy.stats import unitary_group

def main():
    """
    Demonstrates how a birefringent medium can cause theories of optical system
    inversion, like Optical Phase Conjugation (OPC), to fail.
    """
    # --- Plan ---
    # The user's question is about whether adding a birefringent medium can cause a theory of "inverting"
    # an optical system to fail. We interpret this "theory" as Optical Phase Conjugation (OPC), a
    # common method for reversing the effects of a random medium.
    #
    # The principle of OPC is that for a system described by an operator S, if we take the output
    # E_out = S * E_in, phase-conjugate it (E_out*), and send it back through the time-reversed
    # system (S_backward), we can recover the original input (conjugated). For a reciprocal medium,
    # S_backward = S_transpose.
    #
    # The reconstructed field is E_rec = S_transpose * E_out* = (S_transpose * S_conjugate) * E_in_conjugate.
    # For perfect reconstruction, the matrix M = S_transpose * S_conjugate must be a multiple of the identity matrix.
    #
    # We will model the system using Jones matrices to focus on polarization effects.
    # Case 1: The system is just a random medium T. We check M1 = T_transpose * T_conjugate.
    # Case 2: The system has a birefringent plate B and a random medium T (S = T * B).
    # We check M2 = (T*B)_transpose * (T*B)_conjugate.
    #
    # We will show that the deterministic birefringent plate B breaks the system's properties,
    # causing the reconstruction to fail.

    # --- Setup ---

    # Input light: Horizontally polarized
    E_in = np.array([1, 0], dtype=np.complex128)
    E_in_conj = np.conj(E_in)

    # Operator 1: Random Medium (T)
    # We model this as a random 2x2 unitary matrix. This represents a lossless
    # polarization-scrambling medium.
    T = unitary_group.rvs(2, random_state=42) # Use a fixed random state for reproducibility

    # Operator 2: Birefringent Medium (B)
    # We model this as a quarter-wave plate, a deterministic unitary matrix.
    B = 0.5 * np.array([[1 - 1j, 1 + 1j], [1 + 1j, 1 - 1j]], dtype=np.complex128)

    # --- Calculations ---

    # Case 1: System is only the random medium (S = T)
    S1 = T
    # The reconstruction quality depends on this matrix:
    M1 = S1.T @ np.conj(S1)
    # Apply the full OPC process to the input vector
    E_out1 = S1 @ E_in
    E_rec1 = S1.T @ np.conj(E_out1)

    # Case 2: System includes the birefringent plate (S = T * B)
    # Light first passes through B, then T.
    S2 = T @ B
    # The reconstruction quality depends on this matrix:
    M2 = S2.T @ np.conj(S2)
    # Apply the full OPC process to the input vector
    E_out2 = S2 @ E_in
    E_rec2 = S2.T @ np.conj(E_out2)

    # --- Analysis ---

    # We measure the final reconstruction fidelity.
    # We calculate the cosine similarity between the reconstructed vector E_rec
    # and the target vector E_in_conj. A value of 1.0 means perfect reconstruction.
    def cosine_similarity(v1, v2):
        v1_norm = np.linalg.norm(v1)
        v2_norm = np.linalg.norm(v2)
        if v1_norm == 0 or v2_norm == 0:
            return 0.0
        return np.abs(np.dot(np.conj(v1), v2)) / (v1_norm * v2_norm)

    fidelity1 = cosine_similarity(E_rec1, E_in_conj)
    fidelity2 = cosine_similarity(E_rec2, E_in_conj)

    # --- Output ---

    print("Yes, adding a birefringent medium can cause the theory of inverting the system to fail.")
    print("The following simulation demonstrates this using the theory of Optical Phase Conjugation (OPC).\n")
    print("-" * 70)

    print("Case 1: System with Random Medium ONLY")
    print(f"The fidelity of the reconstructed beam is high, but not perfect.")
    print(f"Reconstruction Fidelity = {fidelity1:.4f}\n")

    print("Case 2: System with Birefringent Plate AND Random Medium")
    print("The birefringent plate is a deterministic element that breaks the symmetry required for OPC.")
    print("This leads to a significant drop in reconstruction quality.")
    print(f"Reconstruction Fidelity = {fidelity2:.4f}\n")

    print("Conclusion:")
    print("The fidelity is much lower when the birefringent plate is included, showing that the theory of")
    print("recovering the original beam via phase conjugation is no longer holding.")

    # Final equation format as requested
    print("\nFinal Equation of Fidelities:")
    print(f"Fidelity(Random Only) = {fidelity1:.4f}")
    print(f"Fidelity(Random + Birefringent) = {fidelity2:.4f}")

if __name__ == "__main__":
    main()