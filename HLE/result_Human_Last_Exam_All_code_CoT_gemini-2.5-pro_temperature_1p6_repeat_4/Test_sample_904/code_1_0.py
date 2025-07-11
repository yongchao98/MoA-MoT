import sympy as sp

def solve_fluid_equation():
    """
    This function derives the coefficients A(r) and B(r) for the governing
    linear equation of the interfacial shape ξ(r) in a microfluidic system.

    The derivation is based on the linearized Young-Laplace equation in
    cylindrical coordinates.
    """

    # Define symbolic variables for the derivation
    r = sp.Symbol('r', positive=True) # radial position
    gamma = sp.Symbol('γ', positive=True) # surface tension
    xi = sp.Function('ξ')(r) # interfacial displacement

    # 1. Start with the Young-Laplace equation for an axisymmetric surface.
    # The pressure due to surface tension (capillary pressure) is given by:
    # ΔP_cap = γ * (κ1 + κ2), where κ1 and κ2 are the principal curvatures.
    # For a surface ξ(r) in cylindrical coordinates:
    # κ1 = d²ξ/dr² / (1 + (dξ/dr)²)^(3/2)  (curvature in the r-z plane)
    # κ2 = (dξ/dr) / (r * (1 + (dξ/dr)²)^(1/2)) (curvature in the azimuthal plane)
    
    # 2. Linearize the equation.
    # For small displacements, the slope dξ/dr is small, so (dξ/dr)² ≈ 0.
    # The linearized curvatures are:
    # κ1_lin ≈ d²ξ/dr²
    # κ2_lin ≈ (1/r) * dξ/dr
    # So, the linearized capillary pressure is:
    # ΔP_cap ≈ γ * (d²ξ/dr² + (1/r) * dξ/dr)
    
    # 3. The full equation balances the capillary pressure with the electrostatic pressure,
    # which we can represent as the C(r, ξ) term in the final equation.
    # The problem asks for the form: A(r)*ξ'' + B(r)*ξ' + C(r, ξ) = 0.
    # Our derived equation is: γ*ξ'' + (γ/r)*ξ' - ΔP_electrostatic = 0.
    
    # 4. Identify coefficients A(r) and B(r) by comparing the forms.
    A_r = gamma
    B_r = gamma / r
    
    # --- Output the results ---
    print("Derivation of the governing equation for the fluid interface ξ(r):")
    print("-" * 60)
    print("The shape of the interface is determined by the balance between surface tension and electrostatic pressure.")
    print("The linearized Young-Laplace equation gives the pressure due to surface tension, ΔP_cap.")
    print("\nIn cylindrical coordinates, for small displacements, this is:")
    print("ΔP_cap = γ * (d²ξ/dr² + (1/r) * dξ/dr)")
    print("\nThe full governing equation is ΔP_cap - ΔP_elec = 0, where ΔP_elec is the electrostatic pressure.")
    print("Letting C(r, ξ) = -ΔP_elec, we get the form:")
    print("γ * d²ξ/dr² + (γ/r) * dξ/dr + C(r, ξ) = 0")
    print("\nComparing this to the general form A(r) * d²ξ/dr² + B(r) * dξ/dr + C(r, ξ) = 0, we can identify the coefficients.")
    
    print("\nFinal Result:")
    print(f"The coefficient A(r) is: {A_r}")
    print(f"The coefficient B(r) is: {B_r}")
    
    # Final formatted equation
    print("\nThe final equation with the identified coefficients is:")
    print(f"({A_r}) * d²ξ/dr² + ({B_r}) * dξ/dr + C(r, ξ) = 0")


solve_fluid_equation()