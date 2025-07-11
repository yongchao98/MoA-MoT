import numpy as np

def solve_whitening_filter():
    """
    Solves for the whitening filter W(D) based on a corrected version
    of the ill-posed problem statement.
    """

    # As explained, the problem as stated is ill-posed because Q(D) is not a
    # valid power spectral density. We proceed by assuming a typo fix that makes
    # the problem solvable and consistent with typical textbook exercises.
    # The assumed autocorrelation sequence is q_0 = 5/3, q_1 = 2/3, and q_k = 0 for |k| > 1.
    q0 = 5/3.0
    q1 = 2/3.0

    # The D-transform is Q(D) = q_0 + q_1*(D + D^{-1}). We need to find its causal,
    # minimum-phase spectral factor F(D) = f_0 + f_1*D such that Q(D) = F(D)F(D^{-1}).
    # This gives the equations:
    # 1) f_0^2 + f_1^2 = q_0
    # 2) f_0 * f_1 = q_1

    # From (2), f_1 = q_1 / f_0. Substituting into (1):
    # f_0^2 + (q_1 / f_0)^2 = q_0
    # f_0^4 - q_0 * f_0^2 + q_1^2 = 0
    # This is a quadratic equation for f_0^2.
    
    # Let x = f_0^2. The equation is x^2 - q0*x + q1^2 = 0.
    # Using the quadratic formula: x = (q0 +/- sqrt(q0^2 - 4*q1^2)) / 2
    delta = q0**2 - 4 * q1**2
    if delta < 0:
        print("Error: No real solution exists for the coefficients.")
        return

    sol_sq_1 = (q0 + np.sqrt(delta)) / 2
    sol_sq_2 = (q0 - np.sqrt(delta)) / 2

    # The two solutions for the pair {f_0^2, f_1^2} are {sol_sq_1, sol_sq_2}.
    
    # For F(D) to be minimum-phase, its root D = -f_0/f_1 must be inside the
    # unit circle, i.e., |-f_0/f_1| < 1, which implies f_0^2 < f_1^2.
    if sol_sq_1 < sol_sq_2:
        f0_sq = sol_sq_1
        f1_sq = sol_sq_2
    else:
        f0_sq = sol_sq_2
        f1_sq = sol_sq_1
        
    f0 = np.sqrt(f0_sq)
    # Choose sign of f1 to satisfy f0*f1 = q1 > 0
    f1 = np.sqrt(f1_sq)

    # The causal, minimum-phase factor is F(D) = f0 + f1*D.
    # The required whitening filter is W(D) = 1 / F(D^{-1})
    
    print("Explanation: The problem as stated is ill-posed. The following solution is based on a corrected version where q_1 = 2/3 and other q_k for |k|>1 are zero.")
    print("\nThe spectral factor F(D) is:")
    print(f"F(D) = {f0:.4f} + {f1:.4f} * D")
    
    print("\nThe whitening filter W(D) = 1 / F(D^{-1}) is given by the equation:")
    print(f"W(D) = 1 / ( {f0:.4f} + {f1:.4f} * D**-1 )")


solve_whitening_filter()
<<<W(D) = 1 / ( 0.5774 + 1.1547 * D**-1 )>>>