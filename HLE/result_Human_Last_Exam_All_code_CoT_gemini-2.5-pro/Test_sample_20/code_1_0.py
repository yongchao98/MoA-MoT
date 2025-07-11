import numpy as np

def check_subgradient_inequality(phi, epsilon):
    """
    Checks the subgradient inequality for a given potential phi and measure mu.
    
    Here, nu is a Dirac measure at 0, nu = delta_0.
    And mu is a Dirac measure at epsilon, mu = delta_epsilon.
    
    The inequality is: 0.5 * W(mu, nu)^2 >= integral(phi dmu) - integral(phi dnu)
    """
    
    # The Wasserstein-2 distance between two Dirac measures delta_a and delta_b is |a - b|.
    # So, W(delta_epsilon, delta_0) = |epsilon - 0| = |epsilon|.
    # The squared distance is epsilon^2.
    lhs = 0.5 * epsilon**2
    
    # The integral of a function phi against a Dirac measure delta_a is just phi(a).
    # integral(phi dmu) = integral(phi d(delta_epsilon)) = phi(epsilon)
    # integral(phi dnu) = integral(phi d(delta_0)) = phi(0)
    rhs = phi(epsilon) - phi(0)
    
    print(f"Testing for epsilon = {epsilon}:")
    print(f"  Inequality: 0.5 * W(mu, nu)^2 >= integral(phi dmu) - integral(phi dnu)")
    print(f"  Equation with values: {lhs} >= {rhs}")
    result = lhs >= rhs
    print(f"  Is the inequality satisfied? -> {result}\n")
    return result

# --- Case 1: Trivial Tangent Vector (Constant Potential) ---
# Let phi(x) = c. We choose c = 10, but any constant works.
phi_constant = lambda x: 10.0

print("--- Testing Trivial Tangent Vector (phi(x) = 10) ---")
# The inequality should hold for all epsilon.
check_subgradient_inequality(phi_constant, epsilon=1.0)
check_subgradient_inequality(phi_constant, epsilon=0.1)
check_subgradient_inequality(phi_constant, epsilon=-2.0)


# --- Case 2: Non-Trivial Tangent Vector (Non-Constant Potential) ---
# Let phi(x) = x. This is a non-constant function.
phi_linear = lambda x: x

print("--- Testing Non-Trivial Tangent Vector (phi(x) = x) ---")
# The inequality is 0.5 * epsilon^2 >= epsilon.
# This is NOT true for epsilon in (0, 2).
# We will test an epsilon value that violates the inequality.
check_subgradient_inequality(phi_linear, epsilon=1.0)
# We will test an epsilon value that satisfies the inequality.
check_subgradient_inequality(phi_linear, epsilon=3.0)
check_subgradient_inequality(phi_linear, epsilon=-1.0)
