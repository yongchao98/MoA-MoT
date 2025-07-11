import math

def solve_inductance_change():
    """
    This function explains the derivation and prints the final expression
    for the change in mutual inductance per unit length.
    """
    
    print("### Derivation of the Change in Mutual Inductance ###\n")

    # Step 1: Calculate M1 (Mutual inductance of bare circuits per unit length)
    print("Step 1: Mutual Inductance of Bare Circuits (m₁)")
    print("--------------------------------------------------")
    print("The mutual inductance per unit length, m₁, between two circuits is defined by m₁ = Φ₂₁ / I₁, where Φ₂₁ is the magnetic flux per unit length through circuit 2 due to the current I₁ in circuit 1.")
    print("Circuit 1 consists of two long wires with currents +I and -I, separated by distance h. This forms a magnetic dipole.")
    print("The magnetic flux per unit length from circuit 1 through circuit 2 (a similar pair of wires at a distance d) can be calculated exactly.")
    print("The resulting mutual inductance per unit length is:")
    print("m₁ = (μ₀ / (2 * π)) * ln(1 - (h/d)²)")
    print("Since d > h, the argument of the logarithm is less than 1, making m₁ a negative value.\n")

    # Step 2: Calculate M2 (Mutual inductance with concentrators per unit length)
    print("Step 2: Mutual Inductance with Concentrators (m₂)")
    print("-------------------------------------------------")
    print("Each circuit is enclosed in a cylindrical shell with inner radius R₁ and outer radius R₂.")
    print("The shell material has radial permeability μᵣ → ∞ and angular permeability μᶿ → 0.")
    print("These ideal, anisotropic properties force the magnetic field to be contained entirely within the shell. The material acts as a perfect magnetic shield for an internal source.")
    print("Therefore, the magnetic field produced by circuit 1 and its shell at the location of circuit 2 is zero.")
    print("This means the magnetic flux through circuit 2 is zero (Φ₂₁ = 0), and thus the mutual inductance is zero.")
    print("m₂ = 0\n")

    # Step 3: Calculate the change in mutual inductance (Δm)
    print("Step 3: Change in Mutual Inductance (Δm)")
    print("-----------------------------------------")
    print("The change in mutual inductance per unit length is Δm = m₂ - m₁.")
    print("Δm = 0 - [ (μ₀ / (2 * π)) * ln(1 - (h/d)²) ]")
    print("Δm = - (μ₀ / (2 * π)) * ln(1 - (h/d)²)")
    print("Using the property of logarithms -ln(x) = ln(1/x), we can write this in a positive form:\n")

    # Final Expression
    print("### Final Expression ###")
    print("The change in the mutual inductance per unit length is:")
    
    # Print the equation with components as requested
    final_equation = "ΔM/L = (μ₀ / (2 * π)) * ln(d² / (d² - h²))"
    print(final_equation)
    print("\nWhere:")
    print("  μ₀ is the permeability of free space (4π * 10⁻⁷ H/m)")
    print("  π is the mathematical constant Pi")
    print("  d is the distance between the centers of the two circuits")
    print("  h is the distance between the wires within each circuit")
    print("  ln is the natural logarithm")
    # Output each "number" in the equation
    print("The number '2' appears in the denominator.")

if __name__ == "__main__":
    solve_inductance_change()