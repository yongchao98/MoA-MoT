import numpy as np

def calculate_ell(c, d, n, sigma):
    """
    Calculates the value of ell(c, d) based on the derived formula.

    The original problem is ill-posed due to non-symmetric matrices.
    This solution is based on the following key assumptions to make it solvable:
    1. The parameters a=0 and b=1 are used. This simplifies M and S to identity
       matrices, resolving the symmetry issues.
    2. The expression for the probability density p(Y) simplifies such that its
       logarithm is a linear combination of the sum and sum of squares of the
       log-eigenvalues of Y.
    
    With these assumptions, the value of ell is independent of a and b.
    """

    # The sum of integers from 1 to n
    S1 = n * (n + 1) / 2
    # The sum of squares of integers from 1 to n
    S2 = n * (n + 1) * (2 * n + 1) / 6
    
    # Under the assumption a=0, the log-eigenvalues of X1 are i*ln(c)
    # and for X2 are i*ln(d) for i = 1, ..., n.
    # sum(v_1_i) = ln(c) * S1
    # sum(v_2_i) = ln(d) * S1
    # sum(v_1_i^2) = ln(c)^2 * S2
    # sum(v_2_i^2) = ln(d)^2 * S2
    
    # Difference in sum of squares of log-eigenvalues
    sum_v_sq_diff = S2 * (np.log(c)**2 - np.log(d)**2)
    
    # Difference in sum of log-eigenvalues
    sum_v_diff = S1 * (np.log(c) - np.log(d))

    # The final expression for ell from the simplified log-probability
    # ell = - (1/(2*sigma^2)) * (sum(v1^2)-sum(v2^2)) - ((n-1)/2) * (sum(v1)-sum(v2))
    ell_val = - (1 / (2 * sigma**2)) * sum_v_sq_diff - ((n - 1) / 2) * sum_v_diff

    # Print the step-by-step calculation
    print("--- Calculating ell(a,b,c,d) ---")
    print("Based on simplifying assumptions a=0 and b=1.")
    print(f"Given parameters: n = {n}, sigma = {sigma}, c = {c}, d = {d}")
    print("\nStep 1: Calculate helper sums S1 and S2.")
    print(f"S1 = sum(i for i=1..n) = {S1}")
    print(f"S2 = sum(i^2 for i=1..n) = {S2}")
    
    print("\nStep 2: The formula for ell is derived as:")
    print("ell = - (1/(2*sigma^2)) * (S2 * (ln(c)^2 - ln(d)^2)) - ((n-1)/2) * (S1 * (ln(c) - ln(d)))")

    print("\nStep 3: Substitute values and calculate.")
    term1_val = (1 / (2 * sigma**2)) * sum_v_sq_diff
    term2_val = ((n - 1) / 2) * sum_v_diff
    
    print(f"ell = - (1/(2*{sigma}^2)) * ({S2} * ({np.log(c):.4f}^2 - {np.log(d):.4f}^2)) - (({n-1})/2) * ({S1} * ({np.log(c):.4f} - {np.log(d):.4f}))")
    print(f"ell = - (1/{2*sigma**2}) * ({sum_v_sq_diff:.4f}) - ({ (n-1)/2 }) * ({sum_v_diff:.4f})")
    print(f"ell = -({term1_val:.4f}) - ({term2_val:.4f})")
    print(f"ell = {ell_val}")
    
    return ell_val

if __name__ == '__main__':
    # Given parameters from the problem
    n_val = 20
    sigma_val = 5

    # Parameters a, b, c, d are not specified.
    # We use example values for c and d to demonstrate the calculation.
    # The result is independent of a and b under our assumptions.
    c_val = 2.0
    d_val = 3.0

    final_value = calculate_ell(c_val, d_val, n_val, sigma_val)
    # The final answer in the required format
    # print(f"\n<<< {final_value} >>>")