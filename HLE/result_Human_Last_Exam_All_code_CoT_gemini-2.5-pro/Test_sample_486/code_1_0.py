import numpy as np
from scipy.integrate import quad

def integrand(x, R):
    """
    This is the integrand for the energy integral I(R) after integrating over y and z.
    The integral is I(R) = integral from -R to R of this function dx.
    """
    # sech(t) = 1/cosh(t)
    sech_val = 1 / np.cosh(x / np.sqrt(2))
    return (np.pi / 2) * (R**2 - x**2) * (sech_val**4)

def calculate_energy_integral(R):
    """
    Calculates the energy integral I(R) = integral_{B_R} |nabla u|^2 dV
    for the 1D solution u(x,y,z) = tanh(x/sqrt(2)).
    """
    # We use quad for numerical integration.
    # The integral is symmetric, so we can integrate from 0 to R and multiply by 2.
    result, error = quad(lambda x: integrand(x, R), -R, R)
    return result

def main():
    """
    Main function to perform the analysis.
    """
    print("This script computes the energy integral I(R) for the 1D solution u=tanh(x/sqrt(2))")
    print("and verifies that the energy grows quadratically with R (i.e., a=2).")
    print("-" * 70)
    print(f"{'R':<10}{'I(R)':<25}{'I(R) / R^2':<20}")
    print("-" * 70)

    # We test for a range of R values
    radii = [10, 20, 50, 100, 200]

    for R in radii:
        # Since sech^4 decays very fast, for large R the limits of integration
        # don't need to go to R. The effective support is small.
        # However, using the actual limits [-R, R] is more accurate.
        integral_value = calculate_energy_integral(R)
        
        if R > 0:
            ratio = integral_value / (R**2)
            # The asymptotic growth equation is I(R) ~ C * R^a
            # We are testing for a=2, so we compute C = I(R)/R^2
            # The theoretical value for C is (2*sqrt(2)*pi)/3 ~= 2.9619
            print(f"{R:<10.1f}{integral_value:<25.10f}{ratio:<20.10f}")

    print("-" * 70)
    print("As R increases, the ratio I(R)/R^2 converges to a positive constant.")
    print("This demonstrates that the integral grows like R^2.")
    print("The largest possible value for 'a' is therefore 2.")
    
if __name__ == "__main__":
    main()
