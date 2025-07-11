import textwrap

def find_optimal_catalyst():
    """
    Identifies and describes a leading dual-function catalyst system for
    olefin polymerization and polyolefin hydrogenolysis based on a review
    of current scientific literature.
    """

    # 1. Define Catalyst Components
    # Based on research into single-site catalysts for polymer upcycling.
    metal_group = 4
    metal_element = "Zirconium (Zr)"
    support_formula = "SiO2"
    support_name = "Silica"
    active_form = "Single-site hydride species (Zr-H)"

    # 2. Explain the Rationale
    print("Analysis of Dual-Function Catalyst Systems")
    print("="*45)

    explanation = f"""
The challenge is to find a single-site catalyst active for both creating
C-C bonds (polymerization) and breaking them (hydrogenolysis). A leading
candidate is a Group {metal_group} metal, specifically {metal_element}, supported on a
solid oxide.

The {metal_element} atoms are grafted onto a {support_name} ({support_formula}) support, creating
isolated, well-defined active sites. The oxygen atoms of the support
surface act as ligands, stabilizing the metal center.

This system is 'switchable':
- For Olefin Polymerization: It is active at lower temperatures in the presence of monomers.
- For Polyolefin Hydrogenolysis: At higher temperatures and under Hydrogen (H2)
  pressure, the same sites form hydride species ({active_form}) that can
  efficiently cleave the C-C bonds of the polymer backbone into valuable
  short-chain alkanes.
"""
    print(textwrap.dedent(explanation))

    # 3. Present the Final Combination "Equation"
    # This fulfills the requirement to "output each number in the final equation".
    print("\nOptimal Catalyst Combination:")
    print("="*45)
    print(f"  Metal: {metal_element} (from Group {metal_group})")
    print(f"  Support/Ligand: {support_name} ({support_formula})")


# Execute the function to print the findings
find_optimal_catalyst()
