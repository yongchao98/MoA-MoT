import math

def calculate_fermion_scattering_cross_section():
    """
    Calculates the total cross section for fermion-fermion scattering in the
    given theory, assuming the high-energy limit (E >> m, M).
    """
    
    # In the high-energy limit, the spin-averaged squared matrix element <|M|^2>
    # for t- and u-channel pseudoscalar exchange simplifies to 3 * g^4.
    
    # Let's define the components of the final cross section formula.
    # The formula is of the form: sigma = (A * g^B) / (C * pi * E^D)
    
    # Numerator coefficients
    numerator_constant = 3
    g_power = 4
    
    # Denominator coefficients
    # This comes from the general formula for the cross section:
    # sigma = (1/S) * (1 / (64 * pi^2 * s)) * <|M|^2> * integral(dOmega)
    # S = 2 (symmetry factor for identical final particles)
    # s = 4*E^2 (center-of-mass energy squared)
    # <|M|^2> = 3*g^4
    # integral(dOmega) = 4*pi
    # sigma = (1/2) * (1 / (64 * pi^2 * 4*E^2)) * (3*g^4) * (4*pi)
    # sigma = (1/2) * (3*g^4 * 4*pi) / (256 * pi^2 * E^2)
    # sigma = (3 * g^4) / (128 * pi * E^2)
    
    denominator_constant = 128
    pi_symbol = 'pi'
    E_power = 2
    
    # Print the final equation for the total cross section
    print("The total cross section for the scattering of two fermions in the high-energy limit is:")
    
    # The final equation is sigma = (3 * g^4) / (128 * pi * E^2)
    # We print each numerical component as requested.
    print(f"sigma = ({numerator_constant} * g**{g_power}) / ({denominator_constant} * {pi_symbol} * E**{E_power})")

calculate_fermion_scattering_cross_section()