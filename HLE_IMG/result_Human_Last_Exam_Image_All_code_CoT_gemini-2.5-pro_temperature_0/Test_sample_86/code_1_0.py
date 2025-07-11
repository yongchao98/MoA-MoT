import sympy
from sympy import symbols, sqrt, Rational, pretty_print

def calculate_energy_spectrum():
    """
    Calculates the energy spectrum of a harmonic oscillator with a quartic perturbation,
    based on the first-order self-energy diagram.
    """
    # Define the symbols for the physical quantities
    n = symbols('n', integer=True, nonnegative=True)
    m, u = symbols('m u', real=True, positive=True)
    w0 = symbols('omega_0', real=True, positive=True)
    hbar = symbols('hbar', real=True, positive=True)

    # The mean-field approximation of the interaction term leads to a shift
    # in the squared frequency. The new squared frequency w1^2 is:
    # w1^2 = w0^2 + (u / (2*m)) * <x^2>_0
    # The ground state expectation value <x^2>_0 = hbar / (2*m*w0)
    
    # Calculate the shift in the squared frequency
    delta_w_sq = (u * hbar) / (4 * m**2 * w0)
    
    # Calculate the new squared frequency
    w1_sq = w0**2 + delta_w_sq
    
    # The new frequency is the square root of the new squared frequency
    w1 = sqrt(w1_sq)
    
    # The energy spectrum for a harmonic oscillator with frequency w1 is En = (n + 1/2)*hbar*w1
    En = (n + Rational(1, 2)) * hbar * w1

    # Print the results step-by-step
    print("The self-energy diagram predicts a shift in the oscillator's frequency.")
    print("The shift in the squared frequency (delta_w_sq) is:")
    pretty_print(delta_w_sq)
    print("\nThis leads to a new squared frequency (w1_sq = w0^2 + delta_w_sq):")
    pretty_print(w1_sq)
    print("\nThe new frequency (w1) is:")
    pretty_print(w1)
    print("\nThe predicted energy spectrum En = (n + 1/2) * hbar * w1 is:")
    pretty_print(En)
    
    # For clarity, let's also show the final equation with all numbers explicitly.
    print("\nFinal Equation:")
    print(f"En = (n + 1/2) * hbar * sqrt(omega_0**2 + (1/4) * u * hbar / (m**2 * omega_0))")


calculate_energy_spectrum()