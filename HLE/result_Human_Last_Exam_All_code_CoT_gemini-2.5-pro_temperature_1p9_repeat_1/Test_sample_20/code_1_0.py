import numpy as np
from scipy.linalg import sqrtm

def demonstrate_wasserstein_gradient_at_minimum():
    """
    Demonstrates that the Wasserstein gradient of J(mu) = 0.5 * W_2(mu, nu)^2
    is the zero vector at the minimum (mu=nu).

    This is shown for the case of multivariate Gaussian distributions, where
    the optimal transport map has a known analytical form.
    """
    # Let's use a 2-dimensional space for this example.
    d = 2

    # Define the parameters for the target Gaussian distribution nu.
    mean_nu = np.array([5.0, -2.0])
    # A non-diagonal covariance matrix to make it a bit more general.
    cov_nu = np.array([[2.0, 0.8], [0.8, 1.0]])

    # At the minimum of J, mu must be equal to nu.
    # So, we define the parameters for mu to be identical to nu.
    mean_mu = np.array([5.0, -2.0])
    cov_mu = np.array([[2.0, 0.8], [0.8, 1.0]])

    print("--- Wasserstein Gradient Demonstration ---")
    print("Functional: J(mu) = 0.5 * W_2(mu, nu)^2")
    print("The minimum is at mu = nu.")
    print("\nWe check the gradient at this minimum using multivariate Gaussian distributions.")
    print(f"Target distribution nu: N(mean={mean_nu}, cov={cov_nu.flatten()})")
    print(f"Distribution mu at minimum: N(mean={mean_mu}, cov={cov_mu.flatten()})")

    # For Gaussian distributions, the optimal transport map T: mu -> nu is affine:
    # T(x) = A @ x + b
    # where:
    # A = inv(sqrt(cov_mu)) @ sqrt(sqrt(cov_mu) @ cov_nu @ sqrt(cov_mu)) @ inv(sqrt(cov_mu))
    # b = mean_nu - A @ mean_mu
    
    # We expect that when mu=nu, the map T is the identity map, meaning A=I and b=0.

    # Calculate matrix A
    sqrt_cov_mu = sqrtm(cov_mu)
    inv_sqrt_cov_mu = np.linalg.inv(sqrt_cov_mu)
    inner_term = sqrt_cov_mu @ cov_nu @ sqrt_cov_mu
    sqrt_inner_term = sqrtm(inner_term)
    A = inv_sqrt_cov_mu @ sqrt_inner_term @ inv_sqrt_cov_mu

    # Calculate vector b
    b = mean_nu - A @ mean_mu
    
    # The gradient vector field is v(x) = x - T(x) = x - (A@x + b)
    # The equation for the gradient is: v(x) = (I - A)x - b

    # We check the components of this gradient equation.
    # If the gradient is the zero vector, then (I-A) should be a zero matrix
    # and b should be a zero vector.
    I = np.eye(d)
    grad_linear_part = I - A
    grad_constant_part = -b
    
    print("\n--- The Optimal Transport Map T(x) = Ax + b ---")
    print("The map should be the identity T(x) = x, so A should be the identity matrix and b should be the zero vector.")
    print("\nCalculated matrix A:")
    # print each number in the matrix
    for row in A:
        print(" ".join(f"{num:10.8f}" for num in row))
    
    print("\nCalculated vector b:")
    # print each number in the vector
    print(" ".join(f"{num:10.8f}" for num in b))
    
    print("\n--- The Gradient Equation v(x) = (I - A)x - b ---")
    print("The gradient should be the zero vector, so (I-A) and -b should be zero.")
    
    print("\nLinear part of the gradient (I - A):")
    for row in grad_linear_part:
        print(" ".join(f"{num:10.8f}" for num in row))

    print("\nConstant part of the gradient (-b):")
    print(" ".join(f"{num:10.8f}" for num in grad_constant_part))

    print("\n--- Conclusion ---")
    print("As the numerical results show, A is the identity matrix and b is the zero vector (within floating-point precision).")
    print("Therefore, the gradient v(x) = (I - A)x - b = 0, which is the trivial tangent vector.")
    print("This confirms the theoretical result.")


if __name__ == '__main__':
    demonstrate_wasserstein_gradient_at_minimum()