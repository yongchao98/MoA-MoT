import numpy as np
from scipy.stats import norm

def demonstrate_wasserstein_gradient():
    """
    Demonstrates the concept of the Wasserstein gradient for J(mu) = 0.5 * W(mu, nu)^2.
    The gradient vector field is represented by the displacement map T(x) - x.
    """
    print("--- Wasserstein Gradient Demonstration ---")
    
    # Let's define our target distribution 'nu' as a standard normal distribution.
    # We will use N samples to represent it.
    N = 10
    np.random.seed(0)
    nu_samples = np.sort(norm.rvs(loc=0, scale=1, size=N))
    
    print(f"\nCase 1: mu is NOT equal to nu (J(mu) > 0)")
    # Let 'mu' be a normal distribution shifted by 3.
    mu_samples = np.sort(norm.rvs(loc=3, scale=1, size=N))
    
    # In 1D, the optimal transport map T sends the i-th smallest sample of mu
    # to the i-th smallest sample of nu.
    # T(mu_samples[i]) = nu_samples[i]
    # The gradient is represented by the displacement field: T(x) - x
    # For our discrete samples, this is nu_samples - mu_samples.
    
    displacement_field_1 = nu_samples - mu_samples
    
    print("The distribution 'mu' is represented by samples:")
    print(np.round(mu_samples, 2))
    print("\nThe distribution 'nu' (the target) is represented by samples:")
    print(np.round(nu_samples, 2))
    print("\nThe 'gradient' is the displacement field (nu - mu) that pushes mu to nu:")
    print("This represents the equation for the tangent vector at mu.")
    print(np.round(displacement_field_1, 2))
    print("As you can see, this is a non-trivial (non-zero) tangent vector.")
    
    print("\n" + "="*40 + "\n")
    
    print(f"Case 2: mu IS equal to nu (The minimum of J(mu))")
    # At the minimum of J, we have mu = nu.
    mu_at_minimum_samples = nu_samples.copy()
    
    # Now, let's compute the displacement field for transporting nu to itself.
    # The optimal map T is the identity map, T(x) = x.
    # T(mu_at_minimum_samples[i]) = nu_samples[i]
    
    displacement_field_2 = nu_samples - mu_at_minimum_samples
    
    print("The distribution 'mu' is now at the minimum, so its samples are identical to 'nu':")
    print(np.round(mu_at_minimum_samples, 2))
    print("\nThe 'gradient' is the displacement field (nu - mu) that pushes mu to nu:")
    print("This represents the equation for the tangent vector at the minimum.")
    print(np.round(displacement_field_2, 2))
    print("\nAs predicted, the displacement field is the zero vector.")
    print("This represents the trivial tangent vector at the minimum of J.")

demonstrate_wasserstein_gradient()