import numpy as np

def solve():
    """
    Calculates the trace of the covariance matrix for the given sampling procedure.
    """
    # --- Step 1: Define parameters ---
    d = 101
    alpha = 3.0
    beta = 2.0
    
    # Define vectors v1 and v2 to understand the reflection
    v1 = np.zeros(d)
    v1[0] = 1.0  # e_1
    v2 = np.ones(d)  # 1_d

    # --- Step 2: Calculate E[||v||^2] ---
    # The vector d is constructed as d = [d_1, d_rest], where
    # d_1 = (a-b)/(a+b)
    # d_rest = (2*sqrt(ab) / (||c||*(a+b))) * c
    # The squared norm of d is:
    # ||d||^2 = d_1^2 + ||d_rest||^2
    #         = ((a-b)/(a+b))^2 + (4*a*b)/(a+b)^2
    #         = (a^2 - 2ab + b^2 + 4ab) / (a+b)^2
    #         = (a^2 + 2ab + b^2) / (a+b)^2
    #         = (a+b)^2 / (a+b)^2 = 1.
    # So, d is always a unit vector.
    
    # The vector v is a Householder reflection of d. Reflections are orthogonal
    # transformations, so they preserve the norm.
    # ||v||^2 = ||d||^2 = 1.
    # The expectation of a constant is the constant itself.
    E_norm_v_squared = 1.0

    # --- Step 3: Calculate ||E[v]||^2 ---
    # First, we compute E[d].
    # For d_1 = (a-b)/(a+b):
    # Let Y = a / (a+b). Since a ~ Gamma(alpha, theta) and b ~ Gamma(beta, theta),
    # Y follows a Beta distribution, Y ~ Beta(alpha, beta).
    # d_1 can be rewritten as a/(a+b) - b/(a+b) = Y - (1-Y) = 2*Y - 1.
    # The expectation E[Y] for a Beta(alpha, beta) distribution is alpha / (alpha + beta).
    E_Y = alpha / (alpha + beta)
    E_d1 = 2 * E_Y - 1

    # For the other components d_i (i > 1):
    # d_i involves a factor of c_i, a component of c ~ N(0, I_{d-1}).
    # The distribution of c is spherically symmetric, which means the vector c/||c||
    # is uniformly distributed on the unit sphere. The expectation of any component
    # of c/||c|| is 0. Since a and b are independent of c, E[d_i] becomes a product
    # of expectations, one of which is 0.
    # Thus, E[d_i] = 0 for i > 1.
    E_d = np.zeros(d)
    E_d[0] = E_d1
    
    # Now, compute E[v] = H * E[d], where H is the reflection matrix.
    # H = I - 2 * u * u^T / ||u||^2, where u = v1 - v2.
    # E[v] = H * E[d] = E[d] - (2 / ||u||^2) * u * (u^T * E[d]).
    # Let's check the dot product u^T * E[d].
    # u = v1 - v2 = e1 - 1_d = [0, -1, -1, ..., -1]^T
    # E[d] = [E_d1, 0, 0, ..., 0]^T, which is proportional to e1.
    # u^T * E[d] = u^T * (E_d1 * e1) = E_d1 * (u^T * e1) = E_d1 * 0 = 0.
    # Since the dot product is 0, E[d] lies in the reflection hyperplane and is
    # unaffected by the reflection.
    E_v = E_d

    # Finally, compute the squared norm of E[v].
    norm_E_v_squared = np.dot(E_v, E_v)

    # --- Step 4: Final Calculation ---
    # Tr(Cov(v)) = E[||v||^2] - ||E[v]||^2
    trace_cov_v = E_norm_v_squared - norm_E_v_squared

    # --- Output the final equation as requested ---
    print("The trace of the covariance matrix is calculated as E[||v||^2] - ||E[v]||^2.")
    print("Final Equation:")
    print(f"{E_norm_v_squared} - {norm_E_v_squared} = {trace_cov_v}")

solve()