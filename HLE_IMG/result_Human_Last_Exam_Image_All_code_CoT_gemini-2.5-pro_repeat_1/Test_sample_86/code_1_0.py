import math

def print_energy_spectrum_equation():
    """
    Prints the predicted energy spectrum equation based on the self-energy diagram.
    """
    
    # Define the components of the equation as strings for clear printing
    E_n = "E_n"
    hbar = "hbar"
    n = "n"
    one = 1
    two = 2
    omega_0_sq = "omega_0^2"
    u = "u"
    m_sq = "m^2"
    omega_0 = "omega_0"
    four = 4
    
    # Construct the equation string
    # E_n = hbar * (n + 1/2) * sqrt(omega_0^2 + (u*hbar)/(4*m^2*omega_0))
    
    print("The self-energy diagram predicts a new effective oscillator frequency, omega_new.")
    print("The squared value of this new frequency is:")
    print(f"omega_new^2 = {omega_0_sq} + ({u} * {hbar}) / ({four} * {m_sq} * {omega_0})")
    print("\nThe predicted energy spectrum, E_n, is that of a harmonic oscillator with this new frequency:")
    print(f"{E_n} = {hbar} * ({n} + {one}/{two}) * sqrt({omega_0_sq} + ({u} * {hbar}) / ({four} * {m_sq} * {omega_0}))")

print_energy_spectrum_equation()