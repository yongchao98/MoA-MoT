import sympy

def analyze_soliton_stability():
    """
    Analyzes the stability of a 3D soliton using a scaling argument (Derrick's Theorem).

    The Hamiltonian density is H = A*(grad(m))^2 + D*m.curl(m).
    In 3D, under a scaling of size lambda:
    - The exchange energy (A term) scales as lambda^1.
    - The DMI energy (D term) scales as lambda^2.

    For a localized soliton to form, it must have a non-uniform texture,
    so the initial exchange energy E_ex0 is positive.
    The soliton will also adopt a chiral structure that makes the DMI energy
    negative, so the initial DMI energy E_dmi0 is negative.
    """

    # Define symbolic variables
    # lam is the scaling factor (lambda) for the soliton's size. It must be positive.
    lam = sympy.symbols('lambda', positive=True)

    # E_ex_val represents the positive magnitude of the exchange energy for the unscaled (lambda=1) soliton.
    E_ex_val = sympy.symbols('E_ex0', positive=True)

    # abs_E_dmi_val represents the positive magnitude of the DMI energy for the unscaled (lambda=1) soliton.
    abs_E_dmi_val = sympy.symbols('|E_dmi0|', positive=True)

    print("Step 1: Define the total energy E(lambda) based on scaling laws.")
    # The total energy E is the sum of the scaled components.
    # E(lambda) = (E_ex_val * lambda^1) - (|E_dmi_val| * lambda^2)
    # The numbers 1 and 2 are the scaling exponents for the respective terms.
    E = E_ex_val * lam**1 - abs_E_dmi_val * lam**2
    print("The equation for the total energy of the scaled soliton is:")
    print(f"E(lambda) = {E_ex_val} * lambda**1 - {abs_E_dmi_val} * lambda**2")
    print("-" * 50)

    print("Step 2: Find the energy extremum by solving dE/d(lambda) = 0.")
    # Calculate the first derivative of Energy E with respect to lambda.
    dE_dlam = sympy.diff(E, lam)
    # The derivative equation is dE/d(lambda) = 1 * E_ex0 - 2 * |E_dmi0| * lambda
    print("The first derivative of the energy is:")
    print(f"dE/d(lambda) = {dE_dlam}")

    # Solve for lambda where the derivative is zero to find the extremum.
    lam_extremum_sol = sympy.solve(dE_dlam, lam)
    # Check if a valid solution exists (it should be a single positive value)
    if not lam_extremum_sol:
        print("No energy extremum found for a finite, non-zero size.")
        return
        
    lam_extremum = lam_extremum_sol[0]
    print(f"The energy has an extremum at lambda = {lam_extremum}")
    print("-" * 50)

    print("Step 3: Check the stability by calculating the second derivative d2E/d(lambda)2.")
    # Calculate the second derivative of E with respect to lambda.
    d2E_dlam2 = sympy.diff(dE_dlam, lam)
    # The second derivative is d2E/d(lambda)2 = -2 * |E_dmi0|
    print("The second derivative of the energy is a constant:")
    print(f"d2E/d(lambda)2 = {d2E_dlam2}")
    
    # Analyze the sign of the second derivative.
    # Since abs_E_dmi_val is defined as a positive symbol, the second derivative is always negative.
    is_stable = d2E_dlam2 > 0
    print("\nFor a stable soliton, the second derivative must be positive (a local minimum).")
    print(f"In this case, the second derivative is {d2E_dlam2}, which is always negative.")
    print("Conclusion: The energy extremum is a maximum, not a minimum.")
    print("\nTherefore, any localized soliton described by this Hamiltonian is unstable.")

if __name__ == '__main__':
    analyze_soliton_stability()