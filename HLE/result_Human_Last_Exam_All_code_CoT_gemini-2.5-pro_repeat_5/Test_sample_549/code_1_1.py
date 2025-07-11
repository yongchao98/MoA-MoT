import sympy

def evaluate_quantum_correction():
    """
    This script symbolically evaluates and displays the quantum correction to conductivity
    for an electron in a 3D bulk semiconductor.
    """
    
    # Set up symbols for the equation. These are treated as mathematical variables.
    # e: elementary charge
    # h_bar: reduced Planck constant
    # L_phi: phase coherence length
    # l: elastic mean free path
    # delta_sigma_sym: the symbol for the conductivity correction itself
    e, h_bar, L_phi, l = sympy.symbols('e h_bar L_phi l', real=True, positive=True)
    delta_sigma_sym = sympy.Symbol('delta_sigma')
    pi = sympy.pi

    # Explain the plan and the underlying physics
    print("This script evaluates the quantum correction to conductivity for an electron in a 3D bulk semiconductor.")
    print("The correction, known as the weak localization effect, arises from the constructive interference of an electron's wave function along time-reversed paths.")
    print("This interference increases the probability of the electron returning to its starting point, which enhances backscattering and reduces the overall conductivity.")
    print("\nThe calculation involves integrating a quantity called the 'Cooperon' over all possible closed-loop momenta 'q'.")
    print("-" * 80)
    
    # Explain the key parts of the calculation
    print("Step 1: The Integral Term from Momentum Space.")
    print("The core of the calculation is an integral over the magnitude of momentum, q.")
    print("This integral is evaluated between two physical cutoffs:")
    print(f" - The lower momentum limit is q_min = 1/L_phi, where L_phi is the phase coherence length. Beyond this length, inelastic scattering destroys the interference.")
    print(f" - The upper momentum limit is q_max = 1/l, where l is the elastic mean free path. Below this length, the motion is ballistic, not diffusive.")
    
    print("\nBy convention, the result of this integral part of the formula is written as:")
    integral_term = 1/L_phi - 1/l
    sympy.pprint(integral_term, use_unicode=True)
    print("-" * 80)

    print("Step 2: The Prefactor from Quantum Transport Theory.")
    print("Quantum transport theory provides the prefactor for the 3D case:")
    prefactor = e**2 / (2 * pi**2 * h_bar)
    sympy.pprint(prefactor, use_unicode=True)
    print("-" * 80)

    print("Step 3: Assembling the Final Equation.")
    print("Multiplying the prefactor by the integral term gives the final expression for the quantum correction to conductivity, delta_sigma.")
    
    # Construct the final equation by combining the parts
    final_expression = prefactor * integral_term
    final_equation = sympy.Eq(delta_sigma_sym, final_expression)

    # Print the final equation in a well-formatted way
    print("\nThe final equation is:")
    sympy.pprint(final_equation, use_unicode=True)
    
    # Explicitly print the components of the final equation as requested
    print("\nWhere each term in the final equation is:")
    print(f"-> {delta_sigma_sym}: The quantum correction to conductivity.")
    print(f"-> {e}: The elementary charge.")
    print(f"-> {pi}: The mathematical constant pi (~3.14159).")
    print(f"-> {h_bar}: The reduced Planck constant (h/2pi).")
    print(f"-> {L_phi}: The phase coherence length.")
    print(f"-> {l}: The elastic mean free path.")

if __name__ == '__main__':
    evaluate_quantum_correction()