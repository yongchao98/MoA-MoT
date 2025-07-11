import sympy

def analyze_soliton_stability():
    """
    Analyzes the stability of a localized soliton in a 3D Hamiltonian
    with Heisenberg exchange and Dzyaloshinskii-Moriya interaction (DMI)
    using a scaling argument (Derrick's Theorem).
    """

    # --- Introduction ---
    print("Analyzing the stability of a localized soliton in the 3D Hamiltonian:")
    print("E = integral[ A*(grad(m))^2 + D*m.curl(m) ] dV")
    print("-" * 60)
    print("We use a scaling argument. Let's consider a localized solution m_0(r).")
    print("We create a scaled version m_lambda(r) = m_0(r/lambda), where lambda is a positive scaling factor.")
    print("lambda > 1 corresponds to expansion, lambda < 1 corresponds to contraction.")
    print("We will calculate the energy E(lambda) of this scaled configuration.")
    print("-" * 60)

    # --- Symbolic Scaling Analysis ---
    # Define symbolic variables
    # lambda_ represents the scaling factor lambda
    # E_ex represents the exchange energy integral for the unscaled solution (lambda=1)
    # E_DMI represents the DMI energy integral for the unscaled solution (lambda=1)
    lambda_ = sympy.Symbol('lambda', positive=True)
    E_ex = sympy.Symbol('E_ex', positive=True) # Exchange energy is always positive for a non-uniform soliton
    E_DMI = sympy.Symbol('E_DMI', real=True)

    print("The energy of the scaled configuration can be shown to be:")
    print("E(lambda) = (lambda) * E_ex + (lambda^2) * E_DMI\n")
    print("Where:")
    print("  E_ex = integral[ A*(grad(m_0))^2 ] dV  (Exchange Energy)")
    print("  E_DMI = integral[ D*m_0.curl(m_0) ] dV   (DMI Energy)")
    print("-" * 60)

    # Define the scaled energy function
    E_lambda = lambda_ * E_ex + lambda_**2 * E_DMI

    # --- Condition for a Stationary Solution (Potential Soliton) ---
    print("For a stable soliton to exist at lambda=1, it must be a stationary point of the energy.")
    print("This means the first derivative of E(lambda) with respect to lambda must be zero at lambda=1.")

    # Calculate the first derivative
    dE_dlambda = sympy.diff(E_lambda, lambda_)
    print(f"\nThe first derivative is: dE/d(lambda) = {dE_dlambda}")

    # Evaluate at lambda=1 and set to zero
    equilibrium_eq = sympy.Eq(dE_dlambda.subs(lambda_, 1), 0)
    print("\nSetting the derivative to zero at lambda=1 gives the equilibrium condition:")
    print(f"  {equilibrium_eq.lhs} = {equilibrium_eq.rhs}")
    print("This is the virial-type relation that any stationary solution must satisfy.")
    print("-" * 60)

    # --- Condition for Stability ---
    print("For the soliton to be stable, this stationary point must be an energy minimum.")
    print("This means the second derivative of E(lambda) must be positive at lambda=1.")

    # Calculate the second derivative
    d2E_dlambda2 = sympy.diff(dE_dlambda, lambda_)
    print(f"\nThe second derivative is: d^2E/d(lambda)^2 = {d2E_dlambda2}")

    stability_condition = sympy.Gt(d2E_dlambda2.subs(lambda_, 1), 0)
    print("\nThe stability condition is:")
    print(f"  {stability_condition.lhs} > {stability_condition.rhs}")
    print("-" * 60)

    # --- The Contradiction ---
    print("Now, let's analyze these two conditions together.")
    print(f"1. Equilibrium Condition: {equilibrium_eq.lhs} = 0")
    print(f"2. Stability Condition:  {stability_condition.lhs} > 0")

    print("\nFrom the stability condition (2), we must have E_DMI > 0.")
    print("Now, let's substitute this into the equilibrium condition (1).")
    # We can solve the equilibrium equation for E_ex
    # E_ex = -2 * E_DMI
    print(f"From (1), we can write E_ex = -2 * E_DMI.")
    print("Since we established that E_DMI must be positive for stability, this implies that E_ex must be negative.")
    print("E_ex = -2 * (a positive number) < 0.")

    print("\nHere is the contradiction:")
    print("  - The scaling analysis requires E_ex to be negative for a stable soliton.")
    print("  - However, the exchange energy term, E_ex = integral[ A*(grad(m))^2 ] dV, is always positive (A>0 and (grad(m))^2 >= 0) for any non-uniform (i.e., soliton-like) configuration.")

    print("-" * 60)
    print("Conclusion:")
    print("The conditions for the existence of a stable, localized soliton cannot be simultaneously satisfied with this Hamiltonian.")
    print("Any localized structure would be unstable to either collapsing to a point or expanding indefinitely.")
    print("\nTherefore, it is not possible to stabilize a localized soliton in 3D with only Heisenberg exchange and Dzyaloshinskii-Moriya interaction.")

if __name__ == '__main__':
    analyze_soliton_stability()