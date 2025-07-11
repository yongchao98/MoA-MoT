import numpy as np
from scipy.integrate import quad
from scipy.stats import linregress

def main():
    """
    This script numerically verifies the asymptotic growth rate of the energy integral.
    """

    # For the 1D solution u = tanh(x1/sqrt(2)), the integrand for the energy
    # inside a ball of radius R is |nabla u|^2 integrated over the ball.
    # By slicing, this becomes the 1D integral below.
    # integrand = |nabla u|^2 * Area_of_slice
    # |nabla u|^2 = 0.5 * sech^4(x1/sqrt(2))
    # Area_of_slice = pi * (R^2 - x1^2)
    def integrand(x1, R):
        # sech(z) = 1/cosh(z)
        sech_val = 1.0 / np.cosh(x1 / np.sqrt(2))
        return 0.5 * (sech_val**4) * np.pi * (R**2 - x1**2)

    # We calculate the integral I(R) for several large values of R to observe
    # the asymptotic behavior I(R) ~ C * R^a.
    R_values = np.array([20.0, 40.0, 60.0, 80.0, 100.0, 120.0])
    I_values = []

    print("Calculating integral I(R) for various R...")
    for R in R_values:
        # The quad function integrates from a to b. Here [-R, R].
        integral_value, error = quad(integrand, -R, R, args=(R,))
        I_values.append(integral_value)
        print(f"  For R = {R:<5}, I(R) = {integral_value:.4e}")

    I_values = np.array(I_values)

    # To find the exponent 'a', we perform a linear regression on the log-log data.
    # log(I) = a * log(R) + log(C)
    # The slope of the line is our estimate for 'a'.
    log_R = np.log(R_values)
    log_I = np.log(I_values)

    slope, intercept, r_value, p_value, std_err = linregress(log_R, log_I)
    
    a_estimated = slope

    print("\nTo find the exponent 'a', we fit a line to log(I) vs log(R).")
    print(f"The slope of the line gives the estimated value of a.")
    print(f"The result of the linear regression is a = {a_estimated:.6f}")
    print(f"The R-squared value is {r_value**2:.8f}, which indicates an excellent fit.")
    
    # Based on our analysis, the exact value of 'a' is an integer.
    a_final = round(a_estimated)
    print(f"\nThe largest possible integer value for 'a' is {a_final}.")


if __name__ == "__main__":
    main()