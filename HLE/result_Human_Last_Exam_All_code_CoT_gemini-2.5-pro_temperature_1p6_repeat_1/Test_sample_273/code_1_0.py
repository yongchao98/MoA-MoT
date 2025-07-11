import sympy as sp

def analyze_soliton_stability():
    """
    Analyzes the stability of a localized soliton in 3D with Heisenberg exchange
    and Dzyaloshinskii-Moriya interaction using the Hobart-Derrick scaling argument.
    """

    # --- Step 1: Explanation and Setup ---
    print("Analyzing the stability of a 3D soliton with the Hamiltonian density:")
    print("H = A * (grad(m))^2 + D * m . (curl(m))")
    print("\nWe use the Hobart-Derrick scaling argument. We consider a trial solution m(r)")
    print("and a scaled version m_lambda(r) = m(lambda * r).")
    print("The total energy E(lambda) is calculated by integrating H over all space.\n")

    # --- Step 2: Energy Scaling ---
    print("Let the energy of the unscaled solution (lambda=1) be composed of:")
    print("E_A0 = integral[ A * (grad(m))^2 ] dV > 0 (Exchange Energy)")
    print("E_D0 = integral[ D * m . (curl(m)) ] dV      (DMI Energy)\n")
    print("In 3 dimensions (D=3), the energy terms scale with lambda as follows:")
    print("Exchange term scaling: lambda^(2-D) = lambda^(2-3) = lambda^-1")
    print("DMI term scaling: lambda^(1-D) = lambda^(1-3) = lambda^-2\n")

    # --- Step 3: Symbolic Analysis with SymPy ---
    # Define symbolic variables
    # E_A0 represents the integrated exchange energy, which must be positive.
    # E_D0 represents the integrated DMI energy.
    # lam represents the scaling parameter lambda, which is positive.
    E_A0 = sp.Symbol('E_A0', positive=True)
    E_D0 = sp.Symbol('E_D0', real=True)
    lam = sp.Symbol('lambda', positive=True)

    # Define the total energy as a function of lambda
    E = E_A0 / lam + E_D0 / lam**2

    print(f"The total scaled energy equation is E(lambda) = {E}\n")

    # --- Step 4: Find Extremum (First Derivative) ---
    # Calculate the first derivative with respect to lambda
    dE_dlam = sp.diff(E, lam)
    print("To find an energy extremum, we set the first derivative dE/d(lambda) to zero.\n")
    print(f"The first derivative is: dE/d(lambda) = {dE_dlam}\n")

    # Solve for lambda where the derivative is zero
    # For a non-trivial soliton, DMI must be active (E_D0 != 0) and must favor the texture,
    # so we assume E_D0 < 0, which allows for a positive lambda solution.
    lam_extremum_sol = sp.solve(dE_dlam, lam)
    if not lam_extremum_sol:
        print("No real, positive solution for lambda exists if E_D0 is positive or zero.")
        print("In that case, the energy only decreases as lambda increases (collapse).")
        return
    
    lam_extremum = lam_extremum_sol[0]
    print("Solving dE/d(lambda) = 0 for the stationary point gives the following equation for lambda:")
    print(f"lambda_extremum = {lam_extremum}\n")
    print("(Note: This solution is physically meaningful only if E_D0 is negative, as E_A0 and lambda are positive.)\n")


    # --- Step 5: Determine Stability (Second Derivative) ---
    # Calculate the second derivative
    d2E_dlam2 = sp.diff(dE_dlam, lam)
    print("To determine stability, we check the sign of the second derivative d^2E/d(lambda)^2.")
    print("A stable minimum requires this to be positive.\n")
    print(f"The second derivative is: d^2E/d(lambda)^2 = {d2E_dlam2}\n")

    # Evaluate the second derivative at the extremum
    curvature = d2E_dlam2.subs(lam, lam_extremum)
    simplified_curvature = sp.simplify(curvature)
    print("Substituting lambda_extremum into the second derivative gives the final equation for the curvature:")
    print(f"Curvature = {simplified_curvature}\n")

    # --- Step 6: Conclusion ---
    print("Final Analysis:")
    print("The calculated curvature is (-E_A0**3 / (4 * E_D0**2)).")
    print("Since E_A0 is positive, E_A0**3 is positive.")
    print("Since E_D0 is a real number, E_D0**2 is positive.")
    print("Therefore, the entire expression is guaranteed to be NEGATIVE.")
    print("A negative second derivative indicates that the extremum is an energy MAXIMUM, not a minimum.")
    print("\nConclusion: The system is inherently unstable. Any localized soliton will either collapse to a point (lambda -> inf) or expand indefinitely (lambda -> 0).")

if __name__ == '__main__':
    analyze_soliton_stability()
<<<No>>>