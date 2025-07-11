import numpy as np
from scipy.optimize import fsolve

def solve_elasticity():
    """
    This function calculates the elasticity of the optimal job search intensity q
    with respect to the probability of unemployment p, evaluated at p=0.2.
    """
    # Define the given probability
    p = 0.2

    # --- Step 1: Find the equilibrium q for p=0.2 ---
    # The system of equations can be reduced to a single transcendental equation for q:
    # 4 * (1-q) * exp(2q) - (q^2 - 6q + 7) = 0
    def equilibrium_q_equation(q):
        # We search for a root in the valid domain for q, which is (0, 1).
        if q <= 0 or q >= 1:
            return 1e9 # Return a large value outside the valid range
        return 4 * (1 - q) * np.exp(2 * q) - (q**2 - 6 * q + 7)

    # Solve for the root q_0
    q_0_initial_guess = 0.3
    q_0 = fsolve(equilibrium_q_equation, q_0_initial_guess)[0]

    # --- Step 2: Calculate equilibrium values for b and c_e ---
    # From q = (ln(20) - ln(b))/2  => b = 20 * exp(-2q)
    b_0 = 20 * np.exp(-2 * q_0)
    # From c_e * (1-q) = b * (1.5 - q)
    c_e0 = b_0 * (1.5 - q_0) / (1 - q_0)

    # --- Step 3: Solve for the derivative dq/dp using implicit differentiation ---
    # Implicit differentiation of the system of equations yields a 2x2 linear system
    # for x = dq/dp and z = dc_e/dp of the form:
    # A1*x + B1*z = 0
    # A2*x + B2*z = C2
    
    # Coefficients for the linear system of derivatives
    A1 = -(c_e0 - 4 * b_0 + 2 * b_0 * q_0)
    B1 = 1 - q_0
    A2 = -3 * p * b_0 + 2 * p * b_0 * q_0
    B2 = 1 - p
    C2 = -((20 - c_e0) + (1 - q_0) * b_0)

    # Solve for x = dq/dp
    # From the first equation: z = -(A1/B1) * x
    # Substitute into the second: x * (A2 - B2*A1/B1) = C2
    dq_dp = C2 / (A2 - B2 * A1 / B1)

    # --- Step 4: Calculate the elasticity ---
    elasticity = (dq_dp) * (p / q_0)

    # --- Step 5: Print the results ---
    print("--- Calculation Steps ---")
    print(f"1. At p = {p}, the optimal search intensity is q = {q_0:.4f}")
    print(f"2. The derivative of q with respect to p is dq/dp = {dq_dp:.4f}")
    print("\n--- Final Elasticity Calculation ---")
    print(f"Elasticity E = (dq/dp) * (p / q)")
    print(f"             = {dq_dp:.4f} * ({p} / {q_0:.4f})")
    print(f"             = {elasticity:.3f}")
    
    return elasticity

# Run the calculation and store the final answer
final_answer = solve_elasticity()
print(f"\n<<< {final_answer:.3f} >>>")
