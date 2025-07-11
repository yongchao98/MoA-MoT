import math

def calculate_b22_steric():
    """
    Calculates the second osmotic virial coefficient from steric-only behavior
    for a typical monoclonal antibody (mAb) using the hard-sphere model.

    The steric contribution is an intrinsic property based on the molecule's
    size and molecular weight, and does not depend on solution conditions like
    pH or salt concentration.
    """
    # Standard parameters for a typical monoclonal antibody (mAb)
    # Molecular Weight (Mw) in g/mol
    Mw = 150000.0
    # Steric radius (rs), approximated by hydrodynamic radius, in nm
    rs_nm = 5.0
    # Convert radius from nm to cm for volume calculation (1 nm = 1e-7 cm)
    rs_cm = rs_nm * 1e-7

    # Constants
    # Avogadro's number (Na) in mol^-1
    Na = 6.022e23
    # Pi
    pi = math.pi

    # The formula for the steric-only second osmotic virial coefficient (B22,steric)
    # for a hard sphere is B22,steric = (16 * pi * Na * rs^3) / (3 * Mw)
    # The result will be in cm^3/g, which is equivalent to mL/g.
    b22_steric = (16 * pi * Na * (rs_cm**3)) / (3 * Mw)

    # Print the final equation with all the numbers and the result
    print("The second osmotic virial coefficient from steric-only behavior (B22,steric) is calculated using the hard-sphere model.")
    print("\nAssumed parameters for a typical mAb:")
    print(f"Molecular Weight (Mw): {Mw} g/mol")
    print(f"Steric Radius (rs): {rs_nm} nm")
    print(f"Avogadro's Number (Na): {Na:.3e} /mol")
    
    print("\nFinal Equation:")
    # The f-string formats the numbers for clarity in the equation.
    equation_str = (
        f"B22,steric = (16 * \u03C0 * {Na:.3e} mol\u207B\u00B9 * ({rs_cm:.1e} cm)\u00B3) / "
        f"(3 * {Mw} g/mol)"
    )
    print(equation_str)
    print(f"\nResult: B22,steric = {b22_steric:.3f} mL/g")

calculate_b22_steric()