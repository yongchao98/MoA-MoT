import math

def solve_inductance_change():
    """
    This function explains and derives the expression for the change in mutual inductance
    between two circuits when surrounded by magnetic concentrator shells.
    """
    # Define symbolic representations for the variables
    # Using unicode for better readability in the printout
    mu_0 = "μ₀"
    h = "h"
    pi = "π"
    d = "d"
    R1 = "R₁"
    R2 = "R₂"
    delta_M = "ΔM"
    M1 = "M₁"
    M2 = "M₂"

    print("This script derives the expression for the change in mutual inductance per unit length.")
    print("-" * 70)

    # Step 1: Calculate M1
    print(f"Step 1: Determine the mutual inductance {M1} between the bare circuits.")
    print("In the limit where the distance between circuits 'd' is much larger than the wire separation 'h' (d >> h), each circuit can be modeled as a 2D magnetic dipole.")
    print("The mutual inductance per unit length, M₁, between two such circuits separated by distance 'd' is:")
    print(f"  {M1} = ({mu_0} * {h}²) / (2 * {pi} * {d}²)\n")

    # Step 2: Analyze the concentrator's effect
    print(f"Step 2: Analyze the effect of the concentrator shells on the mutual inductance to find {M2}.")
    print("The circuits are enclosed in shells of a special material with radial permeability μ_r → ∞ and angular permeability μ_θ → 0.")
    print("By solving the magnetostatic boundary value problem for a dipole source inside this type of shell, it's found that the shell amplifies the external magnetic field.")
    print(f"The amplification factor is the ratio of the outer to inner radii: ({R2} / {R1}).")
    print(f"Therefore, the new mutual inductance, {M2}, is the original inductance {M1} amplified by this factor:")
    print(f"  {M2} = {M1} * ({R2} / {R1})\n")

    # Step 3: Calculate the change in inductance
    print(f"Step 3: Calculate the change in mutual inductance ({delta_M}).")
    print(f"The change is defined as the final inductance minus the initial inductance: {delta_M} = {M2} - {M1}.")
    print(f"  {delta_M} = ({M1} * ({R2} / {R1})) - {M1}")
    print(f"  {delta_M} = {M1} * (({R2} / {R1}) - 1)\n")

    # Step 4: Final Expression
    print(f"Step 4: State the final expression for {delta_M}.")
    print(f"Substituting the expression for {M1} from Step 1, we get the final result for the change in mutual inductance per unit length:")
    
    # Print each variable name as requested by the prompt for the final equation.
    print("\nThe final equation is composed of the following variables:")
    print(f"  {mu_0}: Permeability of free space")
    print(f"  {h}: Separation between wires in each circuit")
    print(f"  {d}: Distance between the two circuits")
    print(f"  {R1}: Inner radius of the concentrator shells")
    print(f"  {R2}: Outer radius of the concentrator shells")
    print(f"  {pi}: The mathematical constant pi")
    print("\nFinal Expression:")
    print(f"  {delta_M} = (({mu_0} * {h}²) / (2 * {pi} * {d}²)) * (({R2} / {R1}) - 1)")

solve_inductance_change()