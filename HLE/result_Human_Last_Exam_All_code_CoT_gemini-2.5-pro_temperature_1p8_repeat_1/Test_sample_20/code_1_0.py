import numpy as np

def illustrate_wasserstein_gradient():
    """
    Illustrates that the gradient of the squared Wasserstein distance
    is the trivial (zero) vector at the minimum.
    """
    print("--- Wasserstein Gradient Demonstration ---")
    
    # Define the locations of the target discrete measure 'nu'.
    # For simplicity, we assume uniform weights. Let's use 4 points.
    nu_locs = np.array([-3.0, -1.0, 2.0, 5.0])
    n = len(nu_locs)
    
    print(f"Let nu be a discrete measure with {n} atoms at locations: \n{nu_locs}\n")

    # --- Case 1: mu is not at the minimum (mu != nu) ---
    
    # Define the locations of a different measure 'mu'.
    mu_locs_1 = np.array([-4.0, 0.0, 1.0, 4.5])
    print("Case 1: mu != nu (not at the minimum)")
    print(f"Let mu be a measure with atoms at locations: \n{mu_locs_1}\n")
    
    # In 1D, for sorted locations, the optimal map T takes mu_locs[i] to nu_locs[i].
    # The gradient vector field v(x) is given by x - T(x).
    # Here, we represent it as a vector of displacements for each atom.
    grad_at_mu_1 = mu_locs_1 - nu_locs
    
    print("The gradient vector field at mu is calculated as: Id - T")
    print("On the discrete points, this is a vector v where v_i = mu_locs_1[i] - nu_locs[i].\n")
    print(f"  v = {mu_locs_1} - {nu_locs}")
    print(f"  v = {grad_at_mu_1}\n")
    print("This is a non-zero vector, indicating a direction of descent for J.\n")
    print("-" * 40 + "\n")

    # --- Case 2: mu is at the minimum (mu = nu) ---
    
    mu_locs_2 = nu_locs  # Set mu = nu
    print("Case 2: mu = nu (at the minimum of J)")
    print(f"The measure mu has atoms at the same locations as nu: \n{mu_locs_2}\n")
    
    # The optimal map T from 'nu' to 'nu' is the identity map, T(x) = x.
    # Therefore, T(nu_locs[i]) = nu_locs[i].
    # The gradient vector field is v(x) = x - T(x) = x - x = 0.
    grad_at_mu_2 = mu_locs_2 - nu_locs
    
    print("The gradient vector field at the minimum is calculated as: Id - T_identity")
    print("On the discrete points, this is a vector v where v_i = mu_locs_2[i] - nu_locs[i].\n")
    
    # Print the equation with all numbers
    print(f"  v = {mu_locs_2} - {nu_locs}")
    print(f"  v = {grad_at_mu_2}\n")
    
    print("The result is the zero vector. This is the trivial tangent vector.")

if __name__ == '__main__':
    illustrate_wasserstein_gradient()
