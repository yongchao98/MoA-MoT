def solve_mutual_inductance_change():
    """
    This function derives and prints the expression for the change in mutual
    inductance between two circuits when surrounded by magnetic concentrators.
    """

    # Define the variables using Unicode for better readability in the output
    mu_0 = "μ₀"
    pi = "π"
    h = "h"
    d = "d"
    R1 = "R₁"
    R2 = "R₂"
    M1 = "M₁"
    M2 = "M₂"
    delta_M = "ΔM"

    # --- Step-by-step derivation ---
    print("Here is the step-by-step derivation for the change in mutual inductance (per unit length).")
    print("-" * 75)

    # Step 1: Expression for M1
    print(f"Step 1: Find the mutual inductance of the bare circuits, {M1}.")
    print("The mutual inductance per unit length between two parallel wire pairs, in the limit")
    print(f"where their separation {d} is much larger than the wire spacing {h} ({d} >> {h}), is:")
    print(f"\n  {M1} ≈ -({mu_0} * {h}²) / (2 * {pi} * {d}²)\n")
    print("-" * 75)

    # Step 2: Expression for M2
    print(f"Step 2: Find the mutual inductance with the concentrators, {M2}.")
    print("Each anisotropic magnetic concentrator shell magnifies the external dipole field")
    print(f"of the circuit it encloses by a magnification factor k = {R2}/{R1}.")
    print("Since both the 'source' circuit and the 'receiving' circuit are shielded,")
    print("the mutual inductance is magnified by k * k = k².")
    print(f"\n  {M2} = ({R2}/{R1})² * {M1}\n")
    print("-" * 75)

    # Step 3: Calculating the change ΔM
    print(f"Step 3: Calculate the change {delta_M} = {M2} - {M1}.")
    print("The change in mutual inductance is the difference between the new and original values:")
    print(f"\n  {delta_M} = {M2} - {M1}")
    print(f"  {delta_M} = (({R2}/{R1})² * {M1}) - {M1}")
    print(f"  {delta_M} = (({R2}/{R1})² - 1) * {M1}\n")
    print("-" * 75)

    # Step 4: Final expression for ΔM
    print(f"Step 4: Substitute the expression for {M1} to find the final result for {delta_M}.")
    print("Substituting the approximate value of M1 gives the final expression for the change:")
    print(f"\n  {delta_M} ≈ (({R2}/{R1})² - 1) * [-({mu_0} * {h}²) / (2 * {pi} * {d}²)]\n")

    # Final formatted equation
    print("The final expression for the change in mutual inductance is:")
    final_equation = f"  {delta_M} = -({mu_0} * {h}**2 / (2 * {pi} * {d}**2)) * (({R2}/{R1})**2 - 1)"
    print(final_equation)
    print("-" * 75)

if __name__ == '__main__':
    solve_mutual_inductance_change()