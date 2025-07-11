def solve_fluid_interface_equation():
    """
    Derives and displays the coefficients A(r) and B(r) for the governing linear equation
    of a fluid interface in a cylindrical electrostatic system.
    """

    print("Derivation of the governing linear equation for the interfacial shape ξ(r):")
    print("-------------------------------------------------------------------------")
    print("1. The analysis starts with the Young-Laplace equation, which relates the pressure")
    print("   difference across the interface (ΔP) to the surface tension (γ) and the total")
    print("   curvature (κ_total = κ₁ + κ₂):")
    print("   ΔP = γ * κ_total")
    print("\n2. For an axisymmetric surface z = ξ(r), the total curvature in cylindrical")
    print("   coordinates, linearized for small displacements (where dξ/dr << 1), is:")
    print("   κ_total ≈ d²ξ/dr² + (1/r) * dξ/dr")
    print("\n3. Since gravity is negligible, the pressure difference is due to electrostatic")
    print("   forces (ΔP_electrostatic). The equation is:")
    print("   ΔP_electrostatic = γ * (d²ξ/dr² + (1/r) * dξ/dr)")
    print("\n4. Rearranging this into the standard form A(r) * d²ξ/dr² + B(r) * dξ/dr + C(r, ξ) = 0 gives:")
    print("   (γ) * d²ξ/dr² + (γ/r) * dξ/dr - ΔP_electrostatic = 0")
    print("\n5. By comparing our derived equation to the standard form, we can identify the")
    print("   functions A(r) and B(r):")
    print("\nThe coefficient A(r) is:")
    print("A(r) = γ")
    print("\nThe coefficient B(r) is:")
    print("B(r) = γ / r")

solve_fluid_interface_equation()
