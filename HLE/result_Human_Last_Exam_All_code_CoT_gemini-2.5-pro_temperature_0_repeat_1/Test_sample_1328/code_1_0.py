import numpy as np

def solve_for_optimal_b(P, alpha):
    """
    Calculates the optimal feedback adjustment factor 'b'.

    Args:
        P (float): Half of the total power budget (2P).
        alpha (float): The noise correlation coefficient.
    """
    # The optimal 'b' is a root of the cubic equation:
    # a*b^3 + c*b^2 + d*b + e = 0
    # which is derived from the derivative of the objective function.
    # The equation is: 2*b^3 + 3*alpha*b^2 - (2*P-1)*b - 2*P*alpha = 0

    # Coefficients of the cubic polynomial
    a = 2
    c = 3 * alpha
    d = -(2 * P - 1)
    e = -2 * P * alpha
    
    coeffs = [a, c, d, e]

    print("The optimal 'b' is a real root of the following cubic equation:")
    print(f"{coeffs[0]}*b^3 + {coeffs[1]}*b^2 + {coeffs[2]}*b + {coeffs[3]} = 0\n")

    # Find the roots of the cubic equation
    roots = np.roots(coeffs)

    # Filter for real roots that satisfy the power constraint b^2 <= 2P
    # and check the boundaries of the domain.
    domain_boundary = np.sqrt(2 * P)
    candidate_b = [r.real for r in roots if np.isreal(r) and r.real**2 <= 2 * P]
    candidate_b.extend([-domain_boundary, domain_boundary])
    
    # Objective function to maximize (proportional to mutual information)
    # N(b) = (2P - b^2) * (b^2 + 2*b*alpha + 1)
    def objective_function(b, P, alpha):
        if b**2 > 2 * P:
            return -np.inf # Outside of valid domain
        P1 = 2 * P - b**2
        return P1 * (b**2 + 2 * b * alpha + 1)

    # Find the b that maximizes the objective function
    optimal_b = max(candidate_b, key=lambda b: objective_function(b, P, alpha))

    print(f"For P = {P} and alpha = {alpha}:")
    print(f"The optimal feedback adjustment factor b is: {optimal_b}")

# --- Example Usage ---
# You can change these values to see how the optimal 'b' changes.
P_val = 10.0
alpha_val = 0.5

solve_for_optimal_b(P_val, alpha_val)
