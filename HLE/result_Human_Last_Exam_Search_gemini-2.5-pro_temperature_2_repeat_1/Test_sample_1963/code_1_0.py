import numpy as np
from scipy.optimize import minimize

def get_moments_v(params):
    """Calculates moments based on E and M2."""
    E, M2 = params
    # We need moments up to M_14 for K=7
    M = np.zeros(16)
    M[0] = 1.0
    if M2 <= 0: return None
    M[2] = M2

    # Recursion relation for V(x) = x^2 + x^4
    # (8k+4) M_{2k+2} = (2k-1)(2k-2)(2k-3) M_{2k-4} + 4(2k-1)E M_{2k-2} - 8k M_{2k}
    try:
        # t=1 (corresponds to k=1 in my notes but using direct formula)
        # M_4 = (E - 2*M_2) / 3, is incorrect, the formula is:
        # 12 * M_4 = 4E*M_0 - 8*M_2
        M[4] = (4 * E * M[0] - 8 * M[2]) / 12.0

        # k=2: M_6
        M[6] = ((3*2*1)*M[0] + 4*3*E*M[2] - 8*2*M[4]) / (8*2+4)
        # k=3: M_8
        M[8] = ((5*4*3)*M[2] + 4*5*E*M[4] - 8*3*M[6]) / (8*3+4)
        # k=4: M_10
        M[10] = ((7*6*5)*M[4] + 4*7*E*M[6] - 8*4*M[8]) / (8*4+4)
        # k=5: M_12
        M[12] = ((9*8*7)*M[6] + 4*9*E*M[8] - 8*5*M[10]) / (8*5+4)
        # k=6: M_14
        M[14] = ((11*10*9)*M[8] + 4*11*E*M[10] - 8*6*M[12]) / (8*6+4)
    except OverflowError:
        return None
    return M

def create_hankel_matrices_v(M):
    """Creates the two Hankel matrices for K=7."""
    k_half = 4
    H0 = np.zeros((k_half, k_half))
    H1 = np.zeros((k_half, k_half))
    
    for i in range(k_half):
      for j in range(k_half):
        H0[i, j] = M[2*(i+j)]
        H1[i, j] = M[2*(i+j)+2]
        
    return H0, H1

def cost_function_v(params):
    """Cost function to minimize E with penalty for non-PSD matrices."""
    E, M2 = params
    if E <= 0 or M2 <= 0:
        return 1e12

    M = get_moments_v(params)
    if M is None or np.isinf(M).any() or np.isnan(M).any():
        return 1e12

    H0, H1 = create_hankel_matrices_v(M)
    
    try:
        eigvals_H0 = np.linalg.eigvalsh(H0)
        eigvals_H1 = np.linalg.eigvalsh(H1)
    except np.linalg.LinAlgError:
        return 1e12
        
    penalty = 0.0
    # Add penalty for negative eigenvalues
    penalty += np.sum(np.abs(eigvals_H0[eigvals_H0 < 0]))
    penalty += np.sum(np.abs(eigvals_H1[eigvals_H1 < 0]))
    
    # Cost is E plus a large penalty factor
    return E + 10000 * penalty

# Initial guess based on literature values for similar potentials
initial_guess = np.array([1.39, 0.36])

# Run the optimization
result = minimize(cost_function_v, initial_guess, method='Nelder-Mead', options={'xatol': 1e-8, 'fatol': 1e-8, 'disp': False})

E_opt, M2_opt = result.x[0], result.x[1]

# Display the results
print(f"Minimal E: {E_opt:.3f}")
print(f"Minimal <x^2>: {M2_opt:.3f}")

# Final verification of constraints at the found solution
final_params = (E_opt, M2_opt)
final_moments = get_moments_v(final_params)
H0_final, H1_final = create_hankel_matrices_v(final_moments)
eig_H0 = np.linalg.eigvalsh(H0_final)
eig_H1 = np.linalg.eigvalsh(H1_final)

print(f"\nEigenvalues of H0 at solution: {np.round(eig_H0, 6)}")
print(f"Eigenvalues of H1 at solution: {np.round(eig_H1, 6)}")
print(f"Determinant of H0: {np.linalg.det(H0_final):.2e}")
print(f"Determinant of H1: {np.linalg.det(H1_final):.2e}")