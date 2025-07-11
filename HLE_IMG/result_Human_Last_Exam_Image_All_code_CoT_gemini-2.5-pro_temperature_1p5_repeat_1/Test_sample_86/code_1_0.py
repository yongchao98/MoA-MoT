import sympy
from sympy import symbols, sqrt, pprint

def solve_energy_spectrum():
    """
    Calculates the energy spectrum of a harmonic oscillator with a quartic perturbation,
    based on the first-order self-energy (tadpole) diagram.
    """
    # Define the symbolic variables
    m, w0, u, hbar, n = symbols('m omega_0 u hbar n', real=True, positive=True)
    x = symbols('x')

    # Step 1: Formulate the correction to the squared frequency
    # The new squared frequency w1^2 is related to the old one w0^2 by:
    # m * w1**2 = m * w0**2 + (u/2) * <x^2>
    # where <x^2> is the expectation value for the unperturbed ground state.
    
    # <x^2> for the unperturbed harmonic oscillator ground state
    x_sq_expectation = hbar / (2 * m * w0)

    # Calculate the shift in the squared frequency, Delta(w^2) = w1^2 - w0^2
    # m * (w1^2 - w0^2) = (u/2) * (hbar / (2*m*w0))
    # w1^2 - w0^2 = u*hbar / (4*m^2*w0)
    delta_w_sq = (u * hbar) / (4 * m**2 * w0)

    # Step 2: Calculate the new frequency w1
    w1_sq = w0**2 + delta_w_sq
    w1 = sqrt(w1_sq)

    # Step 3: Determine the new energy spectrum En
    # The self-energy diagram renormalizes the frequency. The predicted spectrum
    # is that of a new harmonic oscillator with frequency w1.
    # The problem asks for the spectrum up to an overall constant E0.
    En = hbar * w1 * (n + sympy.Rational(1, 2))

    # Print the derivation and the final result
    print("The tadpole self-energy diagram leads to a renormalization of the oscillator's frequency.")
    print("Step 1: The correction to the squared frequency (Pi) is calculated.")
    print(f"    Pi = u * hbar / (4 * m**2 * omega_0)")
    
    print("\nStep 2: The new squared frequency is omega_1^2 = omega_0^2 + Pi.")
    print("    omega_1^2 =")
    pprint(w1_sq, use_unicode=True)

    print("\nStep 3: The predicted energy spectrum is E_n = hbar * omega_1 * (n + 1/2) (up to a constant).")
    print("The final expression for the energy spectrum is:")
    print("\nE_n =")
    # Pretty print for better readability
    pprint(En, use_unicode=True)

    # For the final answer format, we will explicitly print the equation components
    print("\nFinal equation in component form:")
    print(f"E_n = {hbar} * sqrt({w0**2} + ({u} * {hbar}) / ({4} * {m**2} * {w0})) * ({n} + {1}/{2})")


solve_energy_spectrum()