import math

def solve_limit():
    """
    This script calculates the limit of n*P(n) as n -> infinity
    based on the analytical solution using the Central Limit Theorem.
    """

    print("Step 1: Define the vectors and the sum S.")
    # v_A = (1, 0), v_B = (0.5, sqrt(3)/2), v_C = (-0.5, sqrt(3)/2)
    # Number of vectors of each type is 2k. Total vectors n = 6k.
    # S = (S_x, S_y) = sum_{i=1 to n} epsilon_i * v_i
    # Let S_A, S_B, S_C be the sums of epsilon_i for each vector type.
    # S_x = 1*S_A + 0.5*S_B - 0.5*S_C
    # S_y = 0*S_A + (sqrt(3)/2)*S_B + (sqrt(3)/2)*S_C
    print("S_x = S_A + 0.5*S_B - 0.5*S_C")
    print("S_y = (sqrt(3)/2)*(S_B + S_C)\n")

    print("Step 2: Apply the Central Limit Theorem.")
    # For large k, S_A, S_B, S_C are approximately normal.
    # E[epsilon_i] = 0, Var(epsilon_i) = 1.
    # E[S_A] = E[S_B] = E[S_C] = 0.
    # S_A, S_B, S_C are sums of 2k Rademacher variables.
    var_S_A = "2k"
    var_S_B = "2k"
    var_S_C = "2k"
    print(f"Var(S_A) = {var_S_A}, Var(S_B) = {var_S_B}, Var(S_C) = {var_S_C}\n")

    print("Step 3: Calculate the covariance matrix of S = (S_x, S_y).")
    # Var(S_x) = Var(S_A) + (0.5)^2*Var(S_B) + (-0.5)^2*Var(S_C)
    # Var(S_x) = 2k + 0.25*(2k) + 0.25*(2k) = 2k + 0.5k + 0.5k = 3k
    var_Sx_factor = 1**2 * 2 + 0.5**2 * 2 + (-0.5)**2 * 2
    
    # Var(S_y) = (sqrt(3)/2)^2 * (Var(S_B) + Var(S_C))
    # Var(S_y) = (3/4) * (2k + 2k) = (3/4) * 4k = 3k
    var_Sy_factor = (math.sqrt(3)/2)**2 * (2 + 2)

    # Cov(S_x, S_y) = E[S_x*S_y] = 0 (as shown in the detailed derivation).
    cov_Sxy_factor = 0
    
    print(f"Var(S_x) = {var_Sx_factor:.1f}k")
    print(f"Var(S_y) = {var_Sy_factor:.1f}k")
    print(f"Cov(S_x, S_y) = {cov_Sxy_factor:.1f}k")
    print("The covariance matrix is [[3k, 0], [0, 3k]].\n")

    print("Step 4: Calculate the probability P(n).")
    # For large k, S is distributed as a 2D Gaussian N(0, [[3k, 0],[0, 3k]]).
    # The PDF is f(x,y) = 1/(2*pi*3k) * exp(-(x^2+y^2)/(2*3k)).
    # P(n) is the integral of this PDF over a disk of radius sqrt(2).
    # P(n) = integral from 0 to 2pi d(theta) * integral from 0 to sqrt(2) [ r * f(r*cos(theta), r*sin(theta)) dr ]
    # P(n) = integral from 0 to sqrt(2) [ (2*pi*r) / (6*pi*k) * exp(-r^2/(6k)) dr ]
    # P(n) = 1/(3k) * integral from 0 to sqrt(2) [ r * exp(-r^2/(6k)) dr ]
    # With substitution u = r^2/(6k), the integral becomes:
    # P(n) = integral from 0 to 2/(6k) [ exp(-u) du ] = 1 - exp(-1/(3k))
    radius_sq = 2
    variance_factor = var_Sx_factor # or var_Sy_factor
    exponent_val_numerator = radius_sq / 2
    exponent_val_denominator = variance_factor
    print(f"P(n) is approximately 1 - exp(-({exponent_val_numerator}/{exponent_val_denominator}k)) = 1 - exp(-1/(3k))\n")
    
    print("Step 5: Evaluate the limit of n*P(n) as n -> infinity.")
    n_factor = 6
    prob_denom_factor = 3
    print(f"The expression to evaluate is: lim_{{k->inf}} {n_factor}k * (1 - exp(-1/({prob_denom_factor}k)))")

    # Let x = 1/(3k). As k -> inf, x -> 0.
    # The expression becomes: lim_{{x->0}} (6k) * (1 - exp(-x))
    # Since 3k = 1/x, 6k = 2/x.
    # lim_{{x->0}} (2/x) * (1 - exp(-x))
    # This is of the form 0/0. Using L'Hopital's rule or Taylor expansion:
    # lim_{{x->0}} 2 * (d/dx (1-exp(-x))) / (d/dx x) = lim_{{x->0}} 2 * exp(-x) / 1 = 2
    
    # The numbers in the final equation lim_{x->0} C * (1 - exp(-x))/x are:
    C = n_factor / prob_denom_factor
    
    print("Let x = 1/(3k). The limit becomes:")
    print(f"lim_{{x->0}} {C:.0f} * (1 - exp(-x))/x")
    print("The final numbers in the equation are C=2, and the limit of (1-exp(-x))/x is 1.")

    final_result = C * 1
    print(f"\nThe value of the limit is {final_result:.1f}.")

solve_limit()
<<<2>>>