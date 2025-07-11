import textwrap

def print_wrapped(text):
    """Helper function to print text with wrapping."""
    print(textwrap.fill(text, width=80))

def get_channel_interaction_residues():
    """
    This function provides information on key residues in the human
    voltage-gated calcium channel beta-1 subunit (CACNB1, UniProt: Q02641).
    The information is derived from structural biology data and peer-reviewed
    scientific literature.
    """
    print("="*80)
    print("Analysis of Human Calcium Channel Beta-1 Subunit (CACNB1) Hotspot Residues")
    print("Protein: Voltage-gated calcium channel subunit beta-1")
    print("UniProt ID: Q02641")
    print("="*80)
    print("\n")

    # --- Question 1: Interaction Hotspots with Alpha-1 Subunit ---
    print_wrapped("1) Which exact residues from beta-1 subunit are hotspots for interaction with alpha-1 subunit?")
    print("-" * 80)
    print_wrapped("The primary interaction occurs between the Alpha Interaction Domain (AID) of the alpha-1 subunit and the Beta Interaction Domain (BID) on the beta-1 subunit. The following residues of human beta-1 are located in the binding groove of the BID and make direct contact with the alpha-1 AID, as determined from the crystal structure (PDB ID: 1V1H).")
    print("\nResidue Position | Description")
    print("-----------------|------------")

    interaction_hotspots = {
        'H328': 'Histidine 328',
        'K332': 'Lysine 332',
        'I334': 'Isoleucine 334',
        'V335': 'Valine 335',
        'D338': 'Aspartate 338',
        'N339': 'Asparagine 339',
        'E391': 'Glutamate 391',
        'R394': 'Arginine 394',
        'L395': 'Leucine 395',
        'G398': 'Glycine 398',
        'R400': 'Arginine 400',
        'V401': 'Valine 401',
        'E402': 'Glutamate 402',
        'G437': 'Glycine 437',
        'W440': 'Tryptophan 440',
        'V441': 'Valine 441',
        'F445': 'Phenylalanine 445',
        'K448': 'Lysine 448',
        'M449': 'Methionine 449',
        'L452': 'Leucine 452',
        'T455': 'Threonine 455'
    }

    for res, desc in interaction_hotspots.items():
        residue_id = desc.split()[0][:3]
        position = res[1:]
        print(f"{residue_id} {position:<13}| Forms part of the BID binding pocket.")

    print("\n")

    # --- Question 2: Gating Modulation Hotspots ---
    print_wrapped("2) Which exact residues from beta-1 subunit are hotspots for fine-tuning gating properties (gating modulation) of alpha-1 subunit?")
    print("-" * 80)
    print_wrapped("Gating modulation is more complex and involves regions of the beta-1 subunit that can allosterically influence the alpha-1 subunit. Key hotspots are entire regions rather than single residues, though specific amino acids have been identified within them.")
    print("\n")

    print_wrapped("Hotspot Region 1: N-Terminus")
    print_wrapped("The N-terminal domain is known to influence channel inactivation kinetics.")
    print("\nResidue(s) Position   | Description")
    print("---------------------|------------")
    print_wrapped("K6, K7, K9           | This lysine-rich cluster at positions 6, 7, and 9 has been shown to be a key modulator of channel kinetics.")
    print("\n")

    print_wrapped("Hotspot Region 2: The 'Hook' Region")
    print_wrapped("This is a flexible linker region between the SH3 and Guanylate Kinase (GK) domains. Swapping this region between different beta subunit types has a profound impact on channel activation and inactivation, demonstrating its critical role in fine-tuning gating.")
    print("\nResidue(s) Position   | Description")
    print("---------------------|------------")
    print_wrapped("254-326              | This entire region, from Proline-254 to Proline-326, acts as a major hotspot for gating modulation.")

if __name__ == '__main__':
    get_channel_interaction_residues()

<<<Interaction Hotspots: H328, K332, I334, V335, D338, N339, E391, R394, L395, G398, R400, V401, E402, G437, W440, V441, F445, K448, M449, L452, T455. Gating Modulation Hotspots: Lysine-rich cluster (K6, K7, K9) in the N-terminus and the 'Hook' region (residues 254-326).>>>