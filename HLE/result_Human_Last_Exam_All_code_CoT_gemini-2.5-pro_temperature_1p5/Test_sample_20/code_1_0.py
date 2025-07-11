import numpy as np
import torch

def solve_and_explain():
    """
    This function explains the solution to the theoretical question and provides
    a numerical illustration for the concept.
    """
    # --- Part 1: Theoretical Explanation ---
    print("Question: Is it true that the Wasserstein regular subgradient of J(mu) = 0.5 * W(mu,nu)^2 is the trivial tangent vector at the minimum of J?")
    print("\nAnswer: Yes, the statement is true.")
    print("\nHere is the step-by-step reasoning:\n")
    
    print("--- Step 1: Identify the minimum of the functional J ---")
    print("The functional is J(mu) = 0.5 * W(mu, nu)^2, where W is the Wasserstein-2 distance.")
    print("By definition, the Wasserstein distance W(mu, nu) is always non-negative, and W(mu, nu) = 0 if and only if mu = nu.")
    print("Therefore, the functional J(mu) has a unique minimum at mu = nu, where its value is J(nu) = 0.")

    print("\n--- Step 2: The subgradient inequality at the minimum ---")
    print("A vector field `v` (which represents a tangent vector) is in the subgradient del J(mu) if for any other measure `rho`, the following inequality holds:")
    print("  J(rho) - J(mu) >= <v, T - Id>_L2(mu)")
    print("where T is the optimal transport map from `mu` to `rho`, and Id is the identity map.")
    print("We evaluate this inequality at the minimum, mu = nu:")
    print("  J(rho) - J(nu) >= <v, T - Id>_L2(nu)")
    print("Substituting J(nu)=0 and J(rho) = 0.5 * W(rho, nu)^2, we get:")
    print("  0.5 * W(rho, nu)^2 >= <v, T - Id>_L2(nu)")
    print("By definition, W(rho, nu)^2 = integral(||T(x) - x||^2 dnu(x)). Let's denote w = T - x.")
    print("The inequality in L2-norm notation is: 0.5 * ||w||^2 >= <v, w>.")

    print("\n--- Step 3: Test the inequality with a specific choice ---")
    print("This inequality must hold for any `w` that represents a tangent vector at `nu`. The subgradient `v` must also be a tangent vector at `nu`.")
    print("We can test the inequality by choosing a tangent vector `w` that is proportional to `v` itself: let w = alpha * v, for a scalar alpha > 0.")
    print("(This is a valid choice, as it corresponds to transporting the measure `nu` along the direction `v`.)")
    print("Substituting w = alpha * v into the inequality gives:")
    print("  0.5 * ||alpha * v||^2 >= <v, alpha * v>")
    
    print("\n--- Step 4: The final equation and conclusion ---")
    print("Let I = ||v||^2 be the squared L2-norm of v. The inequality becomes:")
    print("  0.5 * alpha^2 * I >= alpha * I")
    print("Since `v` is in the subgradient, this must hold for any choice of alpha > 0. Let's pick alpha = 1.0.")
    alpha = 1.0
    factor = 0.5 * (alpha**2) - alpha
    print(f"For alpha = {alpha}, the final equation simplifies to (0.5 * {alpha**2} - {alpha}) * I >= 0, which is:")
    print(f"  {factor} * I >= 0")
    print("Since I = ||v||^2 is the integral of a non-negative function, I must be non-negative (I >= 0).")
    print("The inequality -0.5 * I >= 0 can only be satisfied if I = 0.")
    print("I = ||v||^2 = 0 implies that the vector field `v` is the zero vector almost everywhere with respect to the measure `nu`.")
    print("\nConclusion: The only element in the subgradient of J at its minimum is the zero vector field, which is the trivial tangent vector.")

    # --- Part 2: Numerical Illustration ---
    print("\n\n--- Numerical Illustration in Finite Dimensions ---")
    print("We can demonstrate this with a finite-dimensional analogue. We represent measures `mu` and `nu` with clouds of N=5 points, `y` and `x`.")
    print("The functional J becomes J(y) = (1/2N) * sum(||y_i - x_i||^2). We compute its gradient with respect to `y` at the minimum (y=x).")
    
    # Configuration
    d = 2  # Dimension of the space
    N = 5  # Number of particles

    # 1. Define nu as N points {x_i}. These are fixed.
    x = torch.tensor([[0.5, 0.1], [-0.3, 0.4], [1.2, -1.0], [0.0, 0.8], [-0.5, -0.7]], dtype=torch.float64)
    
    # 2. Define mu as N points {y_i}. This is the variable.
    # We evaluate the gradient at the minimum, so we set y = x.
    y = x.clone().detach().requires_grad_(True)

    # 3. Define the analogous functional J(y)
    J_val = (1.0 / (2.0 * N)) * torch.sum((y - x)**2)

    # 4. Compute the gradient of J with respect to y
    J_val.backward()
    grad_y = y.grad

    # 5. Print the results
    print(f"\nTarget points `x` (representing nu):\n{x.detach().numpy()}")
    print(f"\nPoints `y` (representing mu), at the minimum y=x:\n{y.detach().numpy()}")
    print(f"\nValue of the functional J at the minimum: {J_val.item():.4f}")
    print(f"\nGradient of J w.r.t. y at the minimum:\n{grad_y.numpy()}")
    print("\nAs the numerical example shows, the gradient at the minimum is the zero vector, which corresponds")
    print("to the trivial tangent vector in the continuous setting.")
    
if __name__ == '__main__':
    solve_and_explain()
    print("\n<<<Yes>>>")
