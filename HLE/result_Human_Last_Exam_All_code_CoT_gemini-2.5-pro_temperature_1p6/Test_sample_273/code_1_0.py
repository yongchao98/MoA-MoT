import sympy as sp

def analyze_soliton_stability():
    """
    Analyzes the stability of a localized soliton in a 3D Hamiltonian with
    Heisenberg exchange and Dzyaloshinskii-Moriya (DMI) interaction using a
    scaling argument (Derrick's theorem).
    """

    # --- Setup ---
    # Define symbolic variables
    # lambda is the scaling factor for the soliton size.
    lmbda = sp.Symbol('lambda', real=True, positive=True)
    # H_A is the exchange energy of the unscaled solution (lambda=1).
    # For a non-uniform soliton, H_A must be positive.
    H_A = sp.Symbol('H_A', real=True, positive=True)
    # H_D is the DMI energy of the unscaled solution (lambda=1).
    H_D = sp.Symbol('H_D', real=True)

    # --- Scaling Argument ---
    # The Hamiltonian is H = integral [A(grad m)^2 + D m.(grad x m)] dV
    # In 3D, under a scaling transformation r -> r/lambda:
    # Exchange energy term H_A = A * integral((grad m)^2) dV scales as lambda^(3-2) = lambda^1.
    # DMI energy term H_D = D * integral(m.(grad x m)) dV scales as lambda^(3-1) = lambda^2.
    H = H_A * lmbda + H_D * lmbda**2

    # Calculate first and second derivatives of the energy with respect to lambda
    dH_dlmbda = sp.diff(H, lmbda)
    d2H_dlmbda2 = sp.diff(H, lmbda, 2)

    # --- Analysis and Output ---
    print("Analyzing soliton stability using Derrick's theorem for the 3D Hamiltonian:")
    print("H = integral [A(grad m)^2 + D m.(grad x m)] dV\n")
    print(f"The scaled energy is H(lambda) = {H}\n")

    print("--- Condition 1: Static Solution (Energy Extremum) ---")
    print("For a static solution, the first derivative of energy must be zero at lambda=1.")
    print(f"dH/dlambda = {dH_dlmbda}")
    # Substitute lambda=1 to get the virial theorem
    virial_equation = sp.Eq(dH_dlmbda.subs(lmbda, 1), 0)
    print(f"At lambda=1, this gives the virial equation: {virial_equation.lhs} = 0")
    print(f"The final equation for a static solution, showing the coefficients, is:")
    print(f"(1)*H_A + (2)*H_D = 0\n")

    # From the virial theorem, solve for H_D
    # Since H_A > 0 (for a non-uniform soliton), H_D must be negative.
    print("For a non-trivial soliton, the exchange energy H_A is positive.")
    print("Therefore, Condition 1 requires H_D to be negative (H_D = -H_A / 2).\n")

    print("--- Condition 2: Stability (Energy Minimum) ---")
    print("For the solution to be stable, the second derivative of energy must be positive at lambda=1.")
    print(f"d^2H/dlambda^2 = {d2H_dlmbda2}")
    # The second derivative is constant, so its value at lambda=1 is the same.
    stability_inequality_lhs = d2H_dlmbda2.subs(lmbda, 1)
    print(f"The stability condition is: {stability_inequality_lhs} > 0")
    print("The final equation for stability, showing the coefficient, is:")
    print(f"(2)*H_D > 0\n")

    print("This implies that H_D must be positive for the solution to be stable.\n")

    print("--- Conclusion ---")
    print("We have a contradiction:")
    print("Condition 1 (static solution) requires H_D < 0.")
    print("Condition 2 (stability)       requires H_D > 0.")
    print("\nThese two necessary conditions cannot be satisfied simultaneously.")
    print("Therefore, it is impossible to stabilize a localized soliton in 3D with only Heisenberg exchange and DMI.")

if __name__ == '__main__':
    analyze_soliton_stability()
