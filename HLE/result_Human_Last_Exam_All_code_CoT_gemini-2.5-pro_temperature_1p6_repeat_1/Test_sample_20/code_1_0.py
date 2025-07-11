import numpy as np

def illustrate_subgradient_at_minimum():
    """
    This function provides a numerical illustration for the subgradient of
    J(mu) = 0.5 * W(mu, nu)^2 at its minimum mu=nu.

    We consider a simple 1D discrete measure nu and perturb it along a
    tangent vector w to get mu_t. We then compute the directional
    derivative of J at nu in the direction w, which is given by:
    lim_{t->0} (J(mu_t) - J(nu)) / t

    The theory states this limit equals <v, w> for v in the subgradient.
    Since J(mu_t) is proportional to t^2, the limit is 0. This implies
    <v, w> = 0 for any w, which means v must be the zero vector.
    """
    print("--- Verifying the subgradient of J(mu) = 0.5 * W(mu, nu)^2 at its minimum ---")

    # Define the base measure nu = 0.5 * delta_{-1} + 0.5 * delta_{1}
    nu_support = np.array([-1.0, 1.0])
    nu_weights = np.array([0.5, 0.5])
    print(f"Base measure nu has support {nu_support} and weights {nu_weights}.\n")

    # The minimum of J is at mu = nu, where J(nu) = 0.
    J_nu = 0.0

    # Let's define a small perturbation step 't'
    t = 1e-6

    # --- Test Case 1: Perturbation by a constant vector field w_1 ---
    # w_1(x) = c. This corresponds to a simple translation.
    c = 2.0
    w1_values = np.array([c, c])
    print(f"Test 1: Perturbing along tangent vector w1(x)={c}, represented by values {w1_values}")

    # Perturbed measure mu_t has support x_i + t * w(x_i)
    mu_t1_support = nu_support + t * w1_values
    
    # For 1D measures with identical weights, the squared Wasserstein-2 distance
    # is the weighted sum of squared distances between sorted support points.
    sq_w2_dist1 = np.sum(nu_weights * (mu_t1_support - nu_support)**2)
    J_mu_t1 = 0.5 * sq_w2_dist1
    directional_derivative1 = (J_mu_t1 - J_nu) / t
    
    print(f"For t = {t}:")
    print(f"  - Support of perturbed measure mu_t1: [{mu_t1_support[0]:.7f}, {mu_t1_support[1]:.7f}]")
    print(f"  - Equation: W(mu_t1, nu)^2 = {nu_weights[0]}*({mu_t1_support[0]:.7f} - {nu_support[0]})^2 + {nu_weights[1]}*({mu_t1_support[1]:.7f} - {nu_support[1]})^2 = {sq_w2_dist1:.2e}")
    print(f"  - J(mu_t1) = 0.5 * {sq_w2_dist1:.2e} = {J_mu_t1:.2e}")
    print(f"  - Directional derivative approx: (J(mu_t1) - J(nu)) / t = ({J_mu_t1:.2e} - {J_nu}) / {t} = {directional_derivative1:.2e}")
    print(f"  This value is very close to 0, which implies <v, w1> = 0.\n")
    
    # --- Test Case 2: Perturbation by a scaling vector field w_2 ---
    # w_2(x) = k*x. This corresponds to scaling the measure outwards.
    k = 3.0
    w2_values = k * nu_support
    print(f"Test 2: Perturbing along tangent vector w2(x)={k}*x, represented by values {w2_values}")
    
    # Perturbed measure mu_t has support x_i + t * w(x_i)
    mu_t2_support = nu_support + t * w2_values

    # Calculate squared Wasserstein-2 distance and J
    sq_w2_dist2 = np.sum(nu_weights * (mu_t2_support - nu_support)**2)
    J_mu_t2 = 0.5 * sq_w2_dist2
    directional_derivative2 = (J_mu_t2 - J_nu) / t

    print(f"For t = {t}:")
    print(f"  - Support of perturbed measure mu_t2: [{mu_t2_support[0]:.7f}, {mu_t2_support[1]:.7f}]")
    print(f"  - Equation: W(mu_t2, nu)^2 = {nu_weights[0]}*({mu_t2_support[0]:.7f} - {nu_support[0]})^2 + {nu_weights[1]}*({mu_t2_support[1]:.7f} - {nu_support[1]})^2 = {sq_w2_dist2:.2e}")
    print(f"  - J(mu_t2) = 0.5 * {sq_w2_dist2:.2e} = {J_mu_t2:.2e}")
    print(f"  - Directional derivative approx: (J(mu_t2) - J(nu)) / t = ({J_mu_t2:.2e} - {J_nu}) / {t} = {directional_derivative2:.2e}")
    print(f"  This value is also very close to 0, which implies <v, w2> = 0.\n")
    
    print("Conclusion: Since the directional derivative is 0 for any direction w, the subgradient v must be the zero vector.")

if __name__ == '__main__':
    illustrate_subgradient_at_minimum()