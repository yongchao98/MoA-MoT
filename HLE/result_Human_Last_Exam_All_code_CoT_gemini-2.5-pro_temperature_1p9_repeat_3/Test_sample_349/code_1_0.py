import numpy as np

def analyze_decomposition(z):
    """
    Analyzes the matrix decomposition for a given z.
    
    A = zB - C
    C = zB - A must be Positive Semidefinite (PSD).
    
    This function takes a specific "worst-case" correlation matrix A,
    constructs the corresponding "nice" matrix B, calculates C,
    and checks if C is PSD by inspecting its eigenvalues.
    """
    
    # A specific correlation matrix which is PSD and has a zero eigenvalue.
    A = np.array([
        [1.0, 0.5, 0.5],
        [0.5, 1.0, -0.5],
        [0.5, -0.5, 1.0]
    ])
    
    # B is a "nice" matrix constructed from A. np.arcsin is element-wise.
    # For diagonal elements, A_ii = 1, np.arcsin(1)=pi/2, so B_ii = (2/pi)*(pi/2) = 1.
    B = (2 / np.pi) * np.arcsin(A)
    
    # C is the residual matrix which must be PSD.
    C = z * B - A
    
    # Check for PSD-ness by computing eigenvalues. All should be >= 0.
    eigenvalues = np.linalg.eigvalsh(C)
    
    print(f"--- Analysis for z = {z:.6f} ---")
    print("Matrix A (Correlation Matrix):")
    print(A)
    print("\nMatrix B ('Nice' Matrix):")
    print(B)
    print("\nMatrix C = z*B - A:")
    print(C)
    print("\nEigenvalues of C:")
    print(eigenvalues)
    
    if np.all(eigenvalues >= -1e-9): # Use tolerance for floating point errors
        print("\nConclusion: C is Positive Semidefinite (all eigenvalues >= 0).")
    else:
        print("\nConclusion: C is NOT Positive Semidefinite (has negative eigenvalues).")

    print("\nVerifying the final equation, C = z*B - A, for the first element:")
    print(f"C[0,0] = z * B[0,0] - A[0,0]")
    print(f"{C[0,0]:.6f} = {z:.6f} * {B[0,0]:.6f} - {A[0,0]:.6f}")

# Test with z = pi/2
z_critical = np.pi / 2
analyze_decomposition(z_critical)

print("\n" + "="*50 + "\n")

# Test with z slightly less than pi/2
z_small = np.pi / 2 - 1e-5
analyze_decomposition(z_small)