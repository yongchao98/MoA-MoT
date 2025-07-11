def solve_superconductor_magnetization():
    """
    Derives and prints the analytical expression for the initial magnetization
    curve of a superconducting strip in a perpendicular magnetic field based on
    the critical-state model.
    """
    print("Derivation of the Initial Magnetization Curve M(H)")
    print("======================================================")
    print("We consider a thin, wide superconducting strip (−a ≤ x ≤ a, −b ≤ y ≤ b, with b << a)")
    print("placed in an external magnetic field H applied along the y-axis.")
    print("The superconductor obeys the critical-state model with a constant critical current density Jc.")
    print("\n--- Step 1: Field Penetration Depth ---")
    print("When the field H is applied, shielding currents with density Jc are induced to oppose the field penetration.")
    print("From Ampere's law (∇ × H = J), for this 1D problem, we have dH/dx = ±Jc.")
    print("To screen a positive field H, Jz = +Jc for x > 0 and Jz = -Jc for x < 0.")
    print("Let's analyze the right side (x > 0). The field penetrates from the surface x=a to a depth x_p.")
    print("The field is zero at the 'flux front', located at x = a - x_p.")
    print("We integrate dH/dx = Jc from the flux front (where H=0) to the surface (where H=H_applied):")
    print("∫[from H=0 to H_applied] dH = ∫[from x=a-x_p to a] Jc dx")
    print("H_applied = Jc * [x] from (a-x_p) to a = Jc * (a - (a - x_p))")
    print("This gives the relation between the applied field H and the penetration depth x_p:")
    print("H = Jc * x_p  or  x_p = H / Jc")
    print("This relation holds until the strip is fully penetrated (x_p = a) at the penetration field H_p = a * Jc.\n")

    print("--- Step 2: Magnetic Moment Calculation ---")
    print("The magnetization is a result of the induced shielding currents.")
    print("The magnetic dipole moment per unit length (m_l) is m_l = ∫ (r × J) dA over the cross-section.")
    print("The y-component of the moment is m_y = -∫ x * Jz(x) dA. With dA = 2b * dx (since b << a):")
    print("m_y = -2b * [ ∫[-a to -a+x_p] x*(-Jc) dx + ∫[a-x_p to a] x*(+Jc) dx ]")
    print("By symmetry, the two integrals are identical. The second integral evaluates to:")
    print("Jc * ∫[a-x_p to a] x dx = Jc * [x²/2]_(a-x_p)^a = (Jc/2) * (a² - (a-x_p)²) = (Jc/2) * (2*a*x_p - x_p²)")
    print("So, the total moment per unit length is m_y = -2b * 2 * [(Jc/2) * (2*a*x_p - x_p²)]")
    print("m_y = -2*b*Jc * (2*a*x_p - x_p²)\n")

    print("--- Step 3: Average Magnetization Calculation ---")
    print("Magnetization M is the magnetic moment per unit volume. The volume per unit length is the cross-sectional area A = (2a)*(2b) = 4ab.")
    print("M = m_y / A = [-2*b*Jc * (2*a*x_p - x_p²)] / (4ab)")
    print("M = - (Jc / (2a)) * (2*a*x_p - x_p²)")
    print("Factoring out x_p gives: M = -Jc * x_p * (1 - x_p / (2a))\n")

    print("--- Step 4: Final Expression M(H) ---")
    print("Finally, we substitute the expression for the penetration depth x_p = H / Jc (from Step 1) into the expression for M (from Step 3).")
    print("M(H) = -Jc * (H / Jc) * (1 - (H / Jc) / (2a))")
    print("Simplifying this equation yields the final analytical expression for the initial magnetization curve:\n")
    
    # Print the final formatted equation
    H_str, a_str, Jc_str = "H", "a", "Jc"
    print("######################################################")
    print(f"##  M({H_str}) = -{H_str} * (1 - {H_str} / (2 * {a_str} * {Jc_str}))")
    print(f"##  (valid for 0 ≤ {H_str} ≤ {a_str}*{Jc_str})")
    print("######################################################\n")
    
    # Print the final equation with each component separated, as requested.
    print("The final equation with each component explicitly shown:")
    print("M(H)", "=", "-", "H", "*", "(", "1", "-", "H", "/", "(", "2", "*", "a", "*", "Jc", ")", ")")


if __name__ == "__main__":
    solve_superconductor_magnetization()
