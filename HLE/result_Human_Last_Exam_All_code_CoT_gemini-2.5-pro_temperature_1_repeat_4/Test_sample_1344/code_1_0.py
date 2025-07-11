import textwrap

def suggest_catalyst():
    """
    This function outlines an optimal catalyst system for both polyolefin
    degradation and olefin polymerization.
    """
    # Define the components of the catalyst system
    metal = "Zirconium (Zr)"
    support = "High-surface-area silica (SiO2)"
    ligand_system = "Tris(siloxy) groups derived from the silica surface, (≡SiO)₃"
    active_species = "(≡SiO)₃ZrH (a single-site zirconium hydride)"

    # Print the detailed explanation
    print("Optimal Catalyst System Proposal:")
    print("-" * 35)
    print(f"1.  Metal:   {metal}")
    print(f"2.  Support: {support}")
    print(f"3.  Ligand:  {ligand_system}")
    print(f"4.  Active Species: {active_species}")
    print("-" * 35)

    explanation_title = "\nRationale for Dual Functionality:\n"
    explanation_text = f"""
The proposed catalyst, a single-site zirconium hydride supported on silica, is uniquely positioned to perform both the synthesis and degradation of polyolefins. Its function can be switched by controlling the reaction conditions.

Function 1: Polyolefin Hydrogenolysis (Degradation)
At high temperatures and in the presence of hydrogen gas (H₂), the {active_species} catalyst efficiently breaks down polyolefin chains.
Mechanism: The process involves a β-alkyl elimination from a zirconacyclopentane intermediate, which cleaves a C-C bond in the polymer backbone. The resulting unsaturated fragments are then immediately hydrogenated by the catalyst, producing short, saturated alkanes (e.g., ethane, propane, butane). The silica support prevents the deactivation of the catalyst under these demanding conditions.

Function 2: Olefin Polymerization (Synthesis)
In the presence of an olefin monomer (like ethylene, C₂H₄) and at lower temperatures, the very same catalyst serves as an effective polymerization initiator.
Mechanism: The Zr-H bond in the {active_species} species is highly reactive. An olefin molecule inserts into the Zr-H bond, forming a Zr-Alkyl bond. This is the first step of chain growth. Subsequent olefin molecules rapidly insert into the new Zr-Alkyl bond, propagating the polymer chain. This is a classic mechanism for single-site Ziegler-Natta type catalysts.
"""
    
    wrapped_explanation = textwrap.dedent(explanation_text)
    print(explanation_title)
    print(wrapped_explanation)

    conclusion_title = "\nConclusion:\n"
    conclusion_text = f"""
This {metal}-based system represents an excellent model for a 'circular economy' catalyst. It can create valuable polymers from simple monomers and, at the end of the plastic's life, break it back down into valuable chemical feedstocks or fuels using hydrogen. The key is the single, well-defined active site ({active_species}) that can mediate both bond-forming and bond-breaking reactions depending on the environment.
"""

    wrapped_conclusion = textwrap.dedent(conclusion_text)
    print(conclusion_title)
    print(wrapped_conclusion)


if __name__ == "__main__":
    suggest_catalyst()
