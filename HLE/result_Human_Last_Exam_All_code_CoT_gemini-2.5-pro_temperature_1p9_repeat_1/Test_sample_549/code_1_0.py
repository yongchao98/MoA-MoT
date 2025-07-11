import sympy

def evaluate_quantum_correction():
    """
    This function symbolically derives the quantum correction to conductivity
    in a 3D disordered system (weak localization).
    """
    # 1. Define physical symbols using sympy
    e, h_bar, D, tau_phi, k, l_e = sympy.symbols('e ħ D τ_φ k l_e', real=True, positive=True)

    # L_phi is the phase coherence length
    L_phi = sympy.Symbol('L_φ', real=True, positive=True)

    print("Step 1: The starting formula for quantum correction to conductivity (Δσ) in 3D is:")
    print("Δσ = - (2 * e² * D / ħ) * ∫ d³k/(2π)³ * [1 / (D*k² + 1/τ_φ)]\n")

    # 2. Integrate over the solid angle in 3D momentum space
    # ∫ d³k f(k) -> ∫ 4πk²dk f(k). So, ∫ d³k/(2π)³ -> (1/(2π²)) ∫ k²dk
    angular_integration_factor = 1 / (2 * sympy.pi**2)
    prefactor = - (2 * e**2 * D / h_bar) * angular_integration_factor

    print("Step 2: After integrating over the solid angle, the formula becomes:")
    print("Δσ = [-e²*D / (π²*ħ)] * ∫ [k² / (D*k² + 1/τ_φ)] dk\n")

    # 3. Perform the definite integral over k
    # The integral runs from k=0 up to a cutoff k_max = 1/l_e
    k_max = 1 / l_e
    integrand = k**2 / (D * k**2 + 1/tau_phi)
    integral_result = sympy.integrate(integrand, (k, 0, k_max))

    print(f"Step 3: We integrate over k from 0 to a physical cutoff k_max = 1/l_e.")
    print("The result of the integral ∫[0 to 1/l_e] [k² / (D*k² + 1/τ_φ)] dk is:")
    # We build the result manually for better presentation, which is equivalent to integral_result
    # result = (1/D) * [k - (1/sqrt(D*tau_phi))*atan(k*sqrt(D*tau_phi))] evaluated at k_max
    integral_expression = (1/D) * (k_max - (1/sympy.sqrt(D*tau_phi))*sympy.atan(k_max*sympy.sqrt(D*tau_phi)))
    print(sympy.pretty(integral_expression.subs(k_max, 1/l_e)), "\n")
    
    # 4. Combine and simplify the full expression
    delta_sigma = prefactor * integral_expression
    
    # Substitute L_phi for sqrt(D*tau_phi) for a standard final form
    delta_sigma_final = delta_sigma.subs(sympy.sqrt(D*tau_phi), L_phi)

    print("Step 4: Combining the prefactor with the integral result gives the full correction:")
    print("Δσ =", sympy.simplify(delta_sigma_final), "\n")
    
    # 5. Analyze the result and provide the main physical conclusion
    print("Step 5: Analyze the expression under the common physical condition where the phase")
    print("coherence length is much larger than the elastic mean free path (L_φ >> l_e).")
    print("In this case, atan(L_φ/l_e) approaches π/2.\n")
    
    final_approx = sympy.simplify(delta_sigma_final.subs(sympy.atan(L_phi/l_e), sympy.pi/2))
    print("The correction Δσ simplifies to: ", final_approx, "\n")

    print("This result contains a temperature-independent term, -(e²/(π²*ħ*l_e)), and a")
    print("temperature-dependent term, which is the primary quantity of interest.\n")
    
    # 6. Final Result
    print("=======================================================================")
    print("Final Evaluation: The temperature-dependent quantum correction to conductivity is:")
    print("\n      e²")
    print("Δσ(T) = ——————————")
    print("    2 * π * ħ * L_φ\n")
    print("Where the symbols are:")
    print("  Δσ(T) : Temperature-dependent change in conductivity")
    print("  e     : Elementary charge")
    print("  π     : The constant Pi")
    print("  ħ     : Reduced Planck constant")
    print("  L_φ   : Phase coherence length, L_φ = sqrt(D*τ_φ)")
    print("=======================================================================")

if __name__ == '__main__':
    evaluate_quantum_correction()
