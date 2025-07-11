import sympy

def analyze_soliton_stability():
    """
    Analyzes the stability of a localized soliton in a 3D system
    with Heisenberg exchange and Dzyaloshinskii-Moriya interaction
    using a scaling argument (Derrick's Theorem).
    """
    # Introduction
    print("This script analyzes the stability of a localized soliton in 3D using a scaling argument.")
    print("The Hamiltonian energy is E = integral[ A*(grad m)^2 + D*m.curl(m) ] dV.")
    print("We test for stability against uniform scaling of the spatial coordinates: r -> lambda * r.")
    print("This corresponds to a scaled solution m_lambda(r) = m_0(r/lambda).")
    print("The total energy E(lambda) becomes a function of the scaling parameter lambda.")
    print("-" * 70)

    # 1. Define symbols and energy expressions
    lam = sympy.Symbol('lambda', positive=True, real=True)
    # E_A1 is the exchange energy of the unscaled solution (lambda=1). It must be positive.
    E_A1 = sympy.Symbol('E_A1', positive=True)
    # E_D1 is the DMI energy of the unscaled solution (lambda=1).
    E_D1 = sympy.Symbol('E_D1', real=True)

    print("In 3 dimensions, the energy terms for the scaled solution scale as follows:")
    # In d-dimensions, the exchange term E_A scales as lambda^(d-2). In d=3, E_A scales as lambda^1.
    E_A_lam = E_A1 * lam
    # The DMI term E_D scales as lambda^(d-1). In d=3, E_D scales as lambda^2.
    E_D_lam = E_D1 * lam**2

    print(f"Exchange Energy Term E_A(lambda) = E_A1 * lambda^1")
    print(f"DMI Energy Term      E_D(lambda) = E_D1 * lambda^2")

    # Total energy as a function of lambda
    E_lam = E_A_lam + E_D_lam
    print(f"\nTotal Energy E(lambda) = {E_lam}")
    print("-" * 70)

    # 2. Condition for a stationary solution (a necessary condition for stability)
    print("For a soliton to be a stable solution, its energy must be at a local minimum.")
    print("A necessary condition is that the first derivative of energy with respect to")
    print("the scaling parameter is zero at lambda = 1 (the unscaled solution).")

    # Calculate the first derivative
    dE_dlam = sympy.diff(E_lam, lam)
    print(f"\nFirst derivative dE/d(lambda) = {dE_dlam}")

    # Evaluate the derivative at lambda = 1
    dE_dlam_at_1 = dE_dlam.subs(lam, 1)
    print(f"At lambda=1, dE/d(lambda) = {dE_dlam_at_1}")

    # The stationary condition equation
    stationary_condition_eq = sympy.Eq(dE_dlam_at_1, 0)
    print(f"\nSo, the stationary condition is: {stationary_condition_eq}")
    print("This is the 'final equation' relating the energy components for a stationary soliton.")
    print("The numbers in this equation are:")
    # The equation is 1*E_A1 + 2*E_D1 = 0
    print(1)
    print(2)
    print(0)
    print("-" * 70)

    # 3. Condition for stability (local minimum)
    print("A further condition for stability is that the energy is a local MINIMUM,")
    print("which means the second derivative of energy must be positive at lambda = 1.")

    # Calculate the second derivative
    d2E_dlam2 = sympy.diff(dE_dlam, lam)
    print(f"\nSecond derivative d^2E/d(lambda)^2 = {d2E_dlam2}")

    # Evaluate the second derivative at lambda = 1
    d2E_dlam2_at_1 = d2E_dlam2.subs(lam, 1)
    print(f"At lambda=1, d^2E/d(lambda)^2 = {d2E_dlam2_at_1}")

    stability_condition = sympy.Gt(d2E_dlam2_at_1, 0)
    print(f"\nSo, the stability condition is: {stability_condition}")
    print("-" * 70)

    # 4. Final Analysis and Conclusion
    print("Let's combine the conditions to check for a contradiction.")
    print(f"1. Stationary Condition: {stationary_condition_eq}  =>  E_A1 = -2 * E_D1")
    print(f"2. Stability Condition: {stability_condition}    =>  E_D1 > 0")
    
    print("\nFor a localized soliton to exist, it must be a non-uniform magnetic texture.")
    print("Therefore, its exchange energy, E_A1 = integral(A*(grad m)^2 dV), must be positive (E_A1 > 0).")
    
    print("\nIf E_A1 > 0, then from the stationary condition (1), we must have -2 * E_D1 > 0, which implies E_D1 < 0.")
    print("However, the stability condition (2) requires E_D1 > 0.")

    print("\nConclusion:")
    print("The two necessary conditions are contradictory (E_D1 must be both negative and positive).")
    print("Therefore, it is impossible to stabilize a localized soliton in the given 3D Hamiltonian.")


if __name__ == '__main__':
    analyze_soliton_stability()