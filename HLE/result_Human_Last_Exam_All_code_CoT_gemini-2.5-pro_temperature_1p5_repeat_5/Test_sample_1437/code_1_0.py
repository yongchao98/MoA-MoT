import numpy as np

def calculate_lindhard_limit():
    """
    Calculates the numerical value of the dimensionless Lindhard function
    in the static (omega=0) and long-wavelength (k->0) limit.
    """
    # The formula for the static dimensionless Lindhard function is:
    # f(x) = 0.5 + (1 - x^2) / (4x) * ln((1+x)/(1-x))
    # where x = k / (2*k_F). We need the limit as x -> 0.

    # To find the limit numerically, we use a very small value for x.
    x = 1e-9

    # The formula can be split into two terms.
    # The first term is a constant.
    term1 = 0.5

    # The second term needs to be evaluated in the limit x->0.
    # Analytically, using the Taylor expansion ln((1+y)/(1-y)) ≈ 2y for small y,
    # the second term (1 - x^2)/(4x) * (2x) simplifies to (1/2) * (1-x^2), which is 0.5 as x->0.
    # We calculate it numerically here.
    log_term = np.log((1 + x) / (1 - x))
    x_factor = (1 - x**2) / (4 * x)
    term2 = x_factor * log_term

    # The final result is the sum of the two terms.
    result = term1 + term2

    print("The dimensionless Lindhard function at zero frequency and momentum is calculated from the limit of:")
    print("f(x) = 0.5 + (1 - x^2)/(4*x) * ln((1+x)/(1-x)) as x -> 0")
    print("\nFor a very small x =", x)
    print("The final value is the sum of two terms:")
    print(f"Term 1 = {term1}")
    print(f"Term 2 = {term2:.9f} (approaches 0.5)")
    
    # Printing the "final equation" as requested
    print("\nFinal equation:")
    print(f"{term1} + {term2:.1f} = {result:.1f}")

    print("\nThe numerical value of the Lindhard polarization function, in its dimensionless form, at k=0 and omega=0 is:")
    print(f"{result:.1f}")


calculate_lindhard_limit()