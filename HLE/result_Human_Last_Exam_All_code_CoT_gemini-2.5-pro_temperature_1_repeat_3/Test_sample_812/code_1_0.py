import sympy as sp

def solve_limit_cumulant():
    """
    Calculates the mean and variance of the random variable X,
    then determines the third cumulant of the limiting distribution
    as described by the Central Limit Theorem.
    """
    
    # Define the symbolic variable for integration
    x = sp.Symbol('x')
    
    # Define the probability density function (PDF) of X_i
    pdf = 3 / (2 * (1 + sp.Abs(x))**4)
    
    # --- Step 1: Calculate the mean (mu) of X_i ---
    # Due to the symmetry of the PDF f(x) = f(-x), the mean is 0.
    # We can formally calculate it.
    mu = sp.integrate(x * pdf, (x, -sp.oo, sp.oo))
    
    # --- Step 2: Calculate the variance (sigma^2) of X_i ---
    # Variance = E[X^2] - (E[X])^2
    e_x_squared = sp.integrate(x**2 * pdf, (x, -sp.oo, sp.oo))
    sigma_squared = e_x_squared - mu**2
    
    # --- Step 3: Apply the Central Limit Theorem ---
    # The variable Y_n converges to a Normal distribution N(0, sigma^2).
    # The parameters of this limiting Normal distribution are:
    limit_dist_mean = 0
    limit_dist_variance = sigma_squared
    
    # The cumulants of a Normal distribution N(mu, sigma^2) are kappa_1 = mu,
    # kappa_2 = sigma^2, and kappa_k = 0 for k >= 3.
    third_cumulant = 0
    
    # --- Step 4: Output the results ---
    print(f"The mean of X_i is mu = {mu}.")
    print(f"The variance of X_i is sigma^2 = {sigma_squared}.")
    print(f"\nBy the Central Limit Theorem, the variable Y_n converges in distribution to a Normal distribution.")
    print(f"The parameters of this limiting distribution are mean = {limit_dist_mean} and variance = {limit_dist_variance}.")
    
    # The final equation is kappa_3 = 0. The numbers are 3 and 0.
    print("\nThe third cumulant (kappa_3) of any Normal distribution is 0.")
    print("\nFinal Equation:")
    print(f"kappa_3 = {third_cumulant}")

solve_limit_cumulant()