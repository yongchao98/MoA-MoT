import sympy as sp

def analyze_soliton_stability():
    """
    Analyzes the stability of a localized soliton in 3D using a scaling argument.
    The Hamiltonian density is H = A*(grad m)^2 + D*m.dot(curl m).
    """

    # Define symbolic variables
    # lam is the scaling factor for the soliton size, lambda > 0
    # E_A is the exchange energy for the unscaled solution (lambda=1)
    # E_D is the DMI energy for the unscaled solution (lambda=1)
    lam = sp.symbols('lambda', positive=True)
    E_A, E_D = sp.symbols('E_A E_D', real=True)

    print("Step 1: Analyze how energy terms scale with size `lambda` in 3D.")
    print("------------------------------------------------------------------")

    # Derivation of scaling:
    # E_A = integral[ A * (grad m)^2 ] dV. In 3D, dV ~ lambda^3 and grad ~ 1/lambda.
    # So, E_A(lambda) scales as (1/lambda^2) * lambda^3 = lambda.
    E_A_scaled = E_A * lam
    print(f"The Heisenberg exchange energy term scales as: E_A(lambda) = E_A * lambda")

    # E_D = integral[ D * m.dot(curl m) ] dV. In 3D, dV ~ lambda^3 and curl ~ 1/lambda.
    # So, E_D(lambda) scales as (1/lambda) * lambda^3 = lambda^2.
    E_D_scaled = E_D * lam**2
    print(f"The Dzyaloshinskii-Moriya energy term scales as: E_D(lambda) = E_D * lambda**2")
    print("\n")

    # Total energy as a function of lambda
    E_total = E_A_scaled + E_D_scaled
    print("Step 2: Express total energy E(lambda) and find the equilibrium condition.")
    print("------------------------------------------------------------------------")
    print(f"The total energy of the scaled soliton is E(lambda) = {E_total}")

    # First derivative to find the stationary point (equilibrium)
    dE_dlam = sp.diff(E_total, lam)
    print(f"The first derivative is dE/dlambda = {dE_dlam}")

    # Equilibrium requires dE/dlambda = 0 at lambda=1
    virial_condition_eq = sp.Eq(dE_dlam.subs(lam, 1), 0)
    print(f"For equilibrium at lambda=1, we must have: {virial_condition_eq}")

    print("\nFor a non-trivial soliton (not uniform), the exchange energy E_A must be positive.")
    print(f"From the equilibrium equation E_A + 2*E_D = 0, this implies E_A = -2*E_D.")
    print("Since E_A > 0, the DMI energy E_D must be negative (E_D < 0).")
    print("\n")

    print("Step 3: Check the stability of the equilibrium point.")
    print("----------------------------------------------------")
    # Second derivative to check for a minimum (stability)
    d2E_dlam2 = sp.diff(E_total, lam, 2)
    stability_condition_val = d2E_dlam2.subs(lam, 1)

    print(f"The second derivative is d^2E/dlambda^2 = {d2E_dlam2}")
    print(f"For stability at lambda=1, we need the second derivative to be positive.")
    print(f"The condition is: {stability_condition_val} > 0")
    print(f"This simplifies to the requirement that E_D > 0.")
    print("\n")

    print("Step 4: The Contradiction.")
    print("----------------------------------------------------")
    print("Equilibrium Condition (Step 2) requires: E_D < 0")
    print("Stability Condition   (Step 3) requires: E_D > 0")
    print("\nThese two necessary conditions are contradictory. It is impossible to satisfy both simultaneously.")
    print("Therefore, no stable localized soliton solution can exist for this Hamiltonian in 3D.")

if __name__ == '__main__':
    analyze_soliton_stability()