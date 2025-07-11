import sympy

def analyze_soliton_stability():
    """
    Analyzes the stability of a 3D localized soliton in a system with
    Heisenberg exchange and Dzyaloshinskii-Moriya interaction using a scaling argument.
    """
    # 1. Define symbolic variables
    # lambda_ is the scaling factor for the soliton's size, must be positive.
    lambda_ = sympy.Symbol('lambda', positive=True)
    # E_ex is the integrated exchange energy for a reference soliton (lambda=1).
    # It's always positive for a non-uniform state.
    E_ex = sympy.Symbol('E_ex', positive=True)
    # E_DMI is the integrated DMI energy for a reference soliton (lambda=1).
    # It must be negative to favor a chiral soliton structure.
    E_DMI = sympy.Symbol('E_DMI', negative=True)

    # 2. Explain the model and scaling
    print("--- Analysis of 3D Soliton Stability ---")
    print("The energy functional is E = integral[A*(grad m)^2 + D*m.(curl m)] dV.")
    print("We test the stability of a localized soliton by scaling its size by a factor lambda.")
    print("\nIn 3D, the energy terms scale as follows:")
    print("Exchange Energy: integral[A*(grad m)^2] dV scales as lambda^(3-2) = lambda^1")
    print("DMI Energy: integral[D*m.(curl m)] dV scales as lambda^(3-1) = lambda^2")
    
    # 3. Define the total energy function
    E_lambda = E_ex * lambda_ + E_DMI * lambda_**2
    print(f"\nThe total energy of the scaled soliton is: E(lambda) = E_ex * lambda + E_DMI * lambda^2")
    print("For a chiral soliton to be considered, we need E_ex > 0 and E_DMI < 0.\n")

    # 4. Find extrema by calculating the first derivative
    dE_dlambda = sympy.diff(E_lambda, lambda_)
    print(f"To find energy extrema, we compute the first derivative: dE/dlambda = {dE_dlambda}")
    
    # Solve for lambda where the derivative is zero
    extrema_sols = sympy.solve(dE_dlambda, lambda_)
    lambda_extremum = extrema_sols[0]
    print(f"Setting the derivative to zero gives an extremum at: lambda = {lambda_extremum}\n")

    # 5. Check stability by calculating the second derivative
    d2E_dlambda2 = sympy.diff(dE_dlambda, lambda_)
    print(f"To check for stability, we compute the second derivative: d^2E/dlambda^2 = {d2E_dlambda2}")
    print("For a stable soliton (an energy minimum), this second derivative must be positive.")
    print(f"However, since E_DMI must be negative, the second derivative d^2E/dlambda^2 = 2*E_DMI is always negative.")
    print("This proves the extremum is an energy MAXIMUM, not a minimum.\n")

    # 6. State the conclusion
    print("--- Conclusion ---")
    print("A soliton in this model sits at a maximum of the energy landscape.")
    print("It is fundamentally unstable and will either collapse to the uniform state (lambda -> 0)")
    print("or expand indefinitely (lambda -> infinity) to lower its energy.")
    print("Therefore, stabilizing localized solitons in 3D is not possible with only these two terms.\n")

    # 7. Provide a numerical example
    print("--- Numerical Example ---")
    E_ex_val = 10.0
    E_DMI_val = -5.0
    print(f"Let's use example values: E_ex = {E_ex_val}, E_DMI = {E_DMI_val}")

    lambda_extremum_val = lambda_extremum.subs({E_ex: E_ex_val, E_DMI: E_DMI_val})
    d2E_dlambda2_val = d2E_dlambda2.subs({E_DMI: E_DMI_val})
    
    print(f"The unstable equilibrium size is lambda = -({E_ex_val}) / (2 * {E_DMI_val}) = {lambda_extremum_val}")
    print(f"The second derivative is 2 * ({E_DMI_val}) = {d2E_dlambda2_val}, which is negative, confirming it's a maximum.")
    
    # Calculate the energy at this maximum
    E_at_extremum = E_lambda.subs({E_ex: E_ex_val, E_DMI: E_DMI_val, lambda_: lambda_extremum_val})
    
    print("\nThe energy of this unstable soliton state is calculated by the equation:")
    # The final equation with numbers, as requested
    print(f"E(lambda={lambda_extremum_val}) = ({E_ex_val}) * ({lambda_extremum_val}) + ({E_DMI_val}) * ({lambda_extremum_val})^2 = {E_at_extremum}")
    print("\nThis energy is higher than the energy of the uniform state (E=0), so the soliton is unstable against collapse.")

if __name__ == '__main__':
    analyze_soliton_stability()