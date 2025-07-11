import numpy as np

def solve_trace_covariance():
    """
    This function calculates the trace of the covariance matrix for the given sampling procedure.
    The solution is derived analytically step-by-step and printed.
    """
    
    # --- Problem Parameters ---
    alpha = 3.0
    beta = 2.0

    print("### Analytical Calculation of Tr(Cov(v)) ###")
    print("-" * 50)

    # --- Step 1: Simplification using Matrix Properties ---
    print("Step 1: Simplify the problem using properties of Trace.")
    print("The vector v is derived from d via a Householder reflection: v = H * d.")
    print("A reflection matrix H is orthogonal (H^T * H = I).")
    print("The covariance of v is Cov(v) = Cov(H * d) = H * Cov(d) * H^T.")
    print("Using the cyclic property of trace, Tr(ABC) = Tr(CAB), we find:")
    print("Tr(Cov(v)) = Tr(H^T * H * Cov(d)) = Tr(I * Cov(d)) = Tr(Cov(d)).")
    print("Therefore, the problem reduces to calculating Tr(Cov(d)).")
    print("-" * 50)

    # --- Step 2: Formula for Tr(Cov(d)) ---
    print("Step 2: Apply the formula Tr(Cov(X)) = E[||X||^2] - ||E[X]||^2.")
    print("We will calculate the two terms, E[||d||^2] and ||E[d]||^2, separately.")
    print("-" * 50)

    # --- Step 3: Calculate E[||d||^2] ---
    print("Step 3: Calculate E[||d||^2].")
    print("The vector d has components d_1 = (a-b)/(a+b) and d_{2:d} = (2*sqrt(ab) / (||c||*(a+b))) * c.")
    print("The squared norm ||d||^2 = d_1^2 + ||d_{2:d}||^2.")
    print("A key algebraic simplification shows that ||d||^2 = ((a-b)/(a+b))^2 + (4ab/(a+b)^2) = 1.")
    print("Since the norm is always 1, its expectation is also 1.")
    E_d_norm_sq = 1.0
    print(f"Result: E[||d||^2] = {E_d_norm_sq}")
    print("-" * 50)

    # --- Step 4: Calculate ||E[d]||^2 ---
    print("Step 4: Calculate ||E[d]||^2 by first finding the expectation vector E[d].")
    
    # --- Step 4a: Calculate E[d_1] ---
    print("\n  Part 4a: Calculate E[d_1].")
    print(f"  For d_1 = (a-b)/(a+b), with a~Gamma({alpha},...) and b~Gamma({beta},...), it can be shown that")
    print(f"  E[d_1] = (alpha - beta) / (alpha + beta).")
    E_d1_num = alpha - beta
    E_d1_den = alpha + beta
    E_d1 = E_d1_num / E_d1_den
    print(f"  E[d_1] = ({alpha} - {beta}) / ({alpha} + {beta}) = {E_d1_num} / {E_d1_den} = {E_d1}")

    # --- Step 4b: Calculate E[d_i] for i > 1 ---
    print("\n  Part 4b: Calculate E[d_i] for i > 1.")
    print("  E[d_i] for i > 1 involves the term E[c_j / ||c||].")
    print("  Since c ~ N(0, I), the vector c/||c|| is uniformly distributed on the unit sphere.")
    print("  By symmetry, E[c_j / ||c||] = 0 for all components j.")
    print("  Therefore, E[d_i] = 0 for all i > 1.")
    
    # --- Step 4c: Calculate ||E[d]||^2 ---
    print("\n  Part 4c: Calculate the squared norm of the expectation vector.")
    print(f"  The expectation vector is E[d] = [{E_d1}, 0, 0, ..., 0]^T.")
    E_d_norm_sq_val = E_d1**2
    print(f"  ||E[d]||^2 = E[d_1]^2 + 0 + ... = ({E_d1})^2 = {E_d_norm_sq_val:.4f}")
    print("-" * 50)

    # --- Step 5: Final Calculation ---
    print("Step 5: Combine the results to get the final answer.")
    print("Tr(Cov(v)) = Tr(Cov(d)) = E[||d||^2] - ||E[d]||^2")
    final_result = E_d_norm_sq - E_d_norm_sq_val
    
    print("\nFinal Equation:")
    # This part outputs each number in the final equation as requested.
    final_E_d_norm_sq = 1
    final_E_d_norm_sq_val = ( (alpha - beta) / (alpha + beta) )**2
    print(f"{final_E_d_norm_sq} - {final_E_d_norm_sq_val:.4f} = {final_result}")
    
    print("\nIn fraction form, the calculation is:")
    final_E_d1_num_int = int(E_d1_num)
    final_E_d1_den_int = int(E_d1_den)
    final_numerator = final_E_d1_den_int**2 - final_E_d1_num_int**2
    final_denominator = final_E_d1_den_int**2
    print(f"1 - ({final_E_d1_num_int}/{final_E_d1_den_int})^2 = 1 - {final_E_d1_num_int**2}/{final_E_d1_den_int**2} = {final_numerator}/{final_denominator}")
    
    print("\n--------------------------------------------------")
    print(f"The final numerical answer for Tr(Cov(v)) is: {final_result}")
    print("--------------------------------------------------")
    
    return final_result

if __name__ == '__main__':
    solve_trace_covariance()