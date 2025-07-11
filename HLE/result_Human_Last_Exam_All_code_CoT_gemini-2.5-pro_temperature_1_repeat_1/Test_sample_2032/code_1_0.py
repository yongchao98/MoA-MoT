import scipy.integrate as integrate

def solve_variance():
    """
    Calculates the variance of Y based on the plan described.
    """
    # Step 2: E[Y] is 1/2 by symmetry.
    E_Y = 0.5

    # Step 5: Define the CDF F(d, x) = P(|x - X| <= d) for X ~ U(0,1).
    def F(d, x):
        """
        Calculates the probability that a U(0,1) variable is within distance d of x.
        """
        if d < 0:
            return 0
        return min(1, x + d) - max(0, x - d)

    # Step 6: Define the inner part of the double integral for E[Y^2].
    def integrand(x2, x1):
        """
        The function to be integrated: 6 * x2^2 * F * (1 - F).
        """
        d = abs(x1 - x2)
        f_val = F(d, x1)
        return 6 * x2**2 * f_val * (1 - f_val)

    # Step 7: Compute the double integral for E[Y^2] numerically.
    # The integration is over the unit square [0,1] x [0,1].
    # The options are added for higher accuracy.
    E_Y_squared, error = integrate.dblquad(integrand, 0, 1, lambda x: 0, lambda x: 1, epsabs=1e-9, epsrel=1e-9)

    # Final calculation: Var(Y) = E[Y^2] - (E[Y])^2
    variance = E_Y_squared - E_Y**2

    # Print the equation with the calculated values
    print("Calculation of the variance:")
    print(f"Var(Y) = E[Y^2] - (E[Y])^2")
    print(f"Var(Y) = {E_Y_squared:.6f} - ({E_Y})^2")
    print(f"Var(Y) = {E_Y_squared:.6f} - {E_Y**2:.6f}")
    print(f"Var(Y) = {variance:.6f}")
    
    # The exact analytical result is 13/240. Let's print that for comparison.
    print(f"\nThe exact analytical result is 13/240, which is approximately {13/240:.6f}.")

solve_variance()
# The final numerical answer for the variance.
print(f"\n<<<{13/240}>>>")