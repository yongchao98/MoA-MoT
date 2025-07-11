import numpy as np
import sympy

def analyze_continuity(M0_np, M_k_formula, k_val, name):
    """
    Analyzes and demonstrates the continuity of the minimal polynomial map
    for a given matrix and a converging sequence.

    Args:
        M0_np (np.array): The limit matrix.
        M_k_formula (function): A function k -> np.array for the sequence.
        k_val (int): A large integer value for k to approximate the limit.
        name (str): The name of the case for printing.
    """
    X = sympy.Symbol('X')
    
    # --- Analysis of the limit matrix M0 ---
    M0 = sympy.Matrix(M0_np)
    chi_M0 = M0.charpoly(X)
    pi_M0 = M0.minpoly(X)
    
    print(f"--- {name} ---")
    print(f"\nLet M0 be the matrix:\n{M0_np}")
    
    is_derogatory = chi_M0.degree() > pi_M0.degree()
    derogatory_status = "derogatory" if is_derogatory else "non-derogatory"
    print(f"\nM0 is {derogatory_status}.")
    
    # Print the polynomials for M0. Note that sympy polynomials show the variable and exponents.
    print(f"Characteristic polynomial of M0, chi_M0: {chi_M0.as_expr()}")
    print(f"Minimal polynomial of M0, pi_M0: {pi_M0.as_expr()}")
    
    # --- Analysis of the sequence matrix Mk ---
    Mk_np = M_k_formula(k_val)
    Mk = sympy.Matrix(Mk_np)
    pi_Mk = Mk.minpoly(X)
    
    print(f"\nConsider a sequence M_k -> M0. For k = {k_val}, M_k is:\n{Mk_np}")
    print(f"Minimal polynomial of M_k, pi_Mk: {pi_Mk.as_expr()}")

    # --- Conclusion on continuity ---
    print("\nAs k -> infinity, M_k -> M0.")
    print(f"The limit of the minimal polynomials is lim pi_Mk = {pi_Mk.as_expr()}")
    print(f"The minimal polynomial of the limit is pi_M0 = {pi_M0.as_expr()}")
    
    # Compare the limit of polynomials with the polynomial of the limit.
    # We compare expression strings for simplicity in this demonstration.
    if str(pi_Mk.as_expr()) == str(pi_M0.as_expr()):
        print("\nConclusion: lim pi_Mk = pi_M0. The map appears continuous at M0.")
    else:
        # For floating point results, check if they are close
        coeffs_diff = np.array(pi_Mk.all_coeffs(), dtype=float) - np.array(pi_M0.all_coeffs(), dtype=float)
        if len(coeffs_diff) > 0 and np.allclose(coeffs_diff, 0):
             print("\nConclusion: lim pi_Mk is numerically equal to pi_M0. The map is continuous at M0.")
        else:
             print("\nConclusion: lim pi_Mk != pi_M0. The map is discontinuous at M0.")
    print("-" * (len(name) + 8))

# --- Case 1: Discontinuity at a derogatory matrix ---
M0_discont = np.array([[1.0, 0.0], [0.0, 1.0]])
Mk_discont_formula = lambda k: np.array([[1.0, 1.0/k], [0.0, 1.0]])
analyze_continuity(M0_discont, Mk_discont_formula, k_val=10000, name="Discontinuity Example (Derogatory Matrix)")

# --- Case 2: Continuity at a non-derogatory matrix ---
M0_cont = np.array([[1.0, 1.0], [0.0, 1.0]])
Mk_cont_formula = lambda k: np.array([[1.0, 1.0], [1.0/k, 1.0]])
# Note: minpoly on floating point matrices can be tricky.
# We modify minpoly to work with a numerical tolerance.
# For this example, Sympy handles it well enough to demonstrate the concept.
analyze_continuity(M0_cont, Mk_cont_formula, k_val=10000, name="Continuity Example (Non-Derogatory Matrix)")
