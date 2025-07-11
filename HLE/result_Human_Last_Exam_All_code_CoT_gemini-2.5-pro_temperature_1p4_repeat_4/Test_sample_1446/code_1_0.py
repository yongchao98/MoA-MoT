def solve_critical_exponent_nu():
    """
    Calculates the critical exponent ν within the mean-field approximation
    for a G-symmetric φ⁴ theory.

    This result is a foundational value in the study of critical phenomena and
    is exact for spatial dimensions d greater than the upper critical
    dimension, d_c = 4.
    """
    
    # In mean-field theory, the exponent ν is derived from the Gaussian
    # part of the Landau-Ginzburg-Wilson Hamiltonian, resulting in a simple
    # rational value.
    numerator = 1
    denominator = 2

    # Calculate the value of ν
    nu_value = numerator / denominator

    # Output the final result as an equation
    print("The precise value of the critical exponent ν in the mean-field approximation of a G₄-theory is given by the equation:")
    print(f"{numerator} / {denominator} = {nu_value}")

solve_critical_exponent_nu()
<<<0.5>>>