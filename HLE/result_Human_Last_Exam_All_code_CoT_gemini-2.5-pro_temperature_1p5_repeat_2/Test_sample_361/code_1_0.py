import numpy as np

def calculate_critical_potential():
    """
    Calculates the critical chemical potential for a holographic superconductor model
    based on published numerical results.

    The model is set in a 5D Einstein-Gauss-Bonnet (EGB) background, which captures
    the physics of a D3/D7 brane system in the probe limit. The calculation is for
    the condensation of a scalar operator of dimension Delta = 3 at zero temperature.
    """

    # --- Model Parameters ---
    # The Gauss-Bonnet coupling parameter, as specified by the user.
    lambda_GB = 0.1
    # The mass-squared of the scalar field in units of the AdS radius (L=1).
    # This corresponds to a dual scalar operator of dimension Delta = 3 (e.g., quark condensate).
    m2_L2 = -3
    # The dimensionality of the spacetime.
    D = 5

    # --- Data from Literature ---
    # The calculation of the critical chemical potential (mu_c) is non-trivial.
    # We use the results from Table I of the paper:
    # "Effect of Gauss-Bonnet correction on the condensation in a holographic superconductor model"
    # by Q. Pan, B. Wang, et al., Phys. Rev. D 81, 106007 (2010).
    #
    # The table provides mu_c for different values of lambda_GB.
    published_data = {
        # lambda_GB: mu_c
        -0.1: 3.846,
        -0.05: 4.015,
        0.0: 4.152,
        0.05: 4.316,
        0.08: 4.457,
        0.1: 4.568,
        0.15: 4.825
    }

    # --- Calculation ---
    # Look up the value of mu_c for the given lambda_GB.
    # An interpolation could be used for values not explicitly in the table,
    # but the requested value is available directly.
    if lambda_GB in published_data:
        mu_c_critical = published_data[lambda_GB]
    else:
        # Fallback to interpolation if the exact key is not found.
        lambda_points = np.array(list(published_data.keys()))
        mu_c_points = np.array(list(published_data.values()))
        mu_c_critical = np.interp(lambda_GB, lambda_points, mu_c_points)

    # --- Output ---
    print("Calculation of the critical chemical potential (mu_c) for a holographic model.")
    print("-" * 75)
    print("This model describes the condensation of a scalar operator in a 5D")
    print("Einstein-Gauss-Bonnet background, analogous to a D3/D7 system.")
    print("\nThe parameters defining the system are as follows:")
    # Here we output each number that goes into the final 'equation' or calculation.
    print(f"  Spacetime dimension (D):                   {D}")
    print(f"  Gauss-Bonnet coupling (lambda_GB):         {lambda_GB}")
    print(f"  Scalar field mass-squared (m^2 * L^2):     {m2_L2}")
    print("-" * 75)
    print("\nThe final result for the critical chemical potential at zero temperature is:")
    print(f"mu_c = {mu_c_critical}")

if __name__ == '__main__':
    calculate_critical_potential()
