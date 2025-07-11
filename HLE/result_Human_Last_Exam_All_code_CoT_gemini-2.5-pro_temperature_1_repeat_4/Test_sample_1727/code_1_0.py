import math

def calculate_amplitude_factor():
    """
    Calculates the stable amplitude factor chi based on the derived formula.

    The factor chi relates the surface topography (e_s) to the basal shear
    stress (S_0) and density contrast (Delta_rho * g) via the equation:
    e_s = chi * (S_0 / (Delta_rho * g))

    The derivation from the Stokes equations under the given boundary conditions
    yields the following formula for chi as a function of the dimensionless
    parameter kh, where k is the wavenumber and h is the plate thickness.
    """

    # The derived equation for the amplitude factor chi
    # Let H = kh
    # chi = 2 * sinh(H) / (1 + cosh(H)^2)

    print("The equation for the stable amplitude factor χ is:")
    print("χ = 2 * sinh(kh) / (1 + (cosh(kh))^2)")
    print("\nThis factor relates the surface topography to the driving basal shear stress.")

    # --- Sample Calculation ---
    # Since k and h are not specified, we use a sample value for the
    # dimensionless product kh to demonstrate the calculation.
    kh = 1.0

    print(f"\nFor a sample calculation, let's assume the dimensionless parameter kh = {kh:.1f}\n")

    # Step 1: Calculate sinh(kh) and cosh(kh)
    sinh_kh = math.sinh(kh)
    cosh_kh = math.cosh(kh)

    # Step 2: Calculate the terms in the equation
    numerator = 2 * sinh_kh
    cosh_kh_sq = cosh_kh**2
    denominator = 1 + cosh_kh_sq

    # Step 3: Calculate the final value of chi
    chi = numerator / denominator

    # Output the step-by-step calculation
    print("Calculation Steps:")
    print(f"1. Let H = kh = {kh:.4f}")
    print(f"2. Calculate the hyperbolic functions:")
    print(f"   sinh(H) = {sinh_kh:.4f}")
    print(f"   cosh(H) = {cosh_kh:.4f}")
    print(f"3. Substitute these values into the equation for χ:")
    print(f"   χ = (2 * sinh(H)) / (1 + cosh(H)^2)")
    print(f"   χ = (2 * {sinh_kh:.4f}) / (1 + ({cosh_kh:.4f})^2)")
    print(f"   χ = {numerator:.4f} / (1 + {cosh_kh_sq:.4f})")
    print(f"   χ = {numerator:.4f} / {denominator:.4f}")
    print(f"   χ = {chi:.4f}")


if __name__ == "__main__":
    calculate_amplitude_factor()
