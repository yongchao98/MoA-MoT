import textwrap

def print_wrapped(text, indent=0):
    """Helper function to print wrapped text with indentation."""
    prefix = ' ' * indent
    # Wrap the text to a reasonable width, subsequent lines are indented
    wrapped_text = textwrap.fill(text, width=80, initial_indent=prefix, subsequent_indent=prefix)
    print(wrapped_text)

def get_calcium_channel_hotspots():
    """
    This script identifies and prints key amino acid residues on the human
    voltage-gated calcium channel beta-1 subunit (CACNB1) that are critical
    for its interaction with and modulation of the alpha-1 subunit (Cav2.3).

    The residue numbering for the GK domain is based on the canonical sequence
    of human CACNB1 (UniProt ID: P54289-2).
    The residue numbering for the N-terminus is based on the well-studied
    beta-1b isoform (UniProt ID: P54289-1).
    """

    # --- 1. Hotspots for Interaction with Alpha-1 Subunit ---
    # These residues are located in the Guanylate Kinase-like (GK) domain of the
    # beta-1 subunit and form a hydrophobic pocket that binds the Alpha Interaction
    # Domain (AID) of the alpha-1 subunit. This binding is the primary structural
    # basis for the alpha-beta subunit assembly.

    interaction_hotspots = {
        "W237 (Tryptophan)": "A key aromatic residue lining the hydrophobic binding groove. Its bulky side chain makes significant contact with the alpha-1 subunit's AID helix.",
        "V296 (Valine)": "Contributes to the hydrophobic surface of the binding pocket, stabilizing the interaction with nonpolar residues on the AID.",
        "I300 (Isoleucine)": "Another crucial hydrophobic residue that deepens the binding pocket and increases the affinity for the AID.",
        "W339 (Tryptophan)": "A highly conserved and critical residue. Mutation of this tryptophan often completely abolishes the high-affinity interaction between the beta and alpha-1 subunits.",
        "V400 (Valine)": "Forms part of the binding interface, contributing to the overall hydrophobic character required for stable complex formation.",
        "L403 (Leucine)": "This residue is located on the surface of the GK domain and participates in the hydrophobic interactions that anchor the alpha-1 subunit."
    }

    # --- 2. Hotspots for Fine-Tuning Gating Properties ---
    # Gating modulation is a complex process. It is initiated by the core interaction
    # (using the residues above), but fine-tuned by other regions of the beta-1 subunit.

    gating_modulation_hotspots = {
        "GK Domain Residues (e.g., W237, W339)": "The fundamental binding of the GK domain to the alpha-1 I-II loop is the primary mechanism for gating modulation. This interaction stabilizes the loop, altering channel activation and inactivation kinetics.",
        "N-terminal K10, K13, K14 (Lysines) in beta-1b isoform": "This polybasic '3K' motif in the N-terminus of the beta-1b isoform (P54289-1) is critical for fine-tuning. It interacts with negatively charged phospholipids in the cell membrane, which allosterically influences channel gating and membrane trafficking."
    }

    print("=" * 80)
    print("      Analysis of Human Calcium Channel Beta-1 Subunit (CACNB1) Hotspots")
    print("=" * 80)
    print("\n")

    # --- Print Interaction Hotspots ---
    print("1) Hotspots for Interaction with Alpha-1 Subunit (Cav2.3)")
    print("-" * 55)
    print_wrapped("The following residues in the Guanylate Kinase-like (GK) domain of the beta-1 subunit are critical for forming the physical binding interface with the alpha-1 subunit:", indent=2)
    print("")
    for residue, description in interaction_hotspots.items():
        print(f"  - Residue: {residue}")
        print_wrapped(f"    Role: {description}", indent=4)
        print("")

    # --- Print Gating Modulation Hotspots ---
    print("\n2) Hotspots for Fine-Tuning Gating Properties of Alpha-1 Subunit")
    print("-" * 68)
    print_wrapped("Gating modulation relies on the core interaction but is fine-tuned by distinct elements. Key hotspots include:", indent=2)
    print("")
    for element, description in gating_modulation_hotspots.items():
        print(f"  - Element: {element}")
        print_wrapped(f"    Role: {description}", indent=4)
        print("")

if __name__ == "__main__":
    get_calcium_channel_hotspots()
    print("<<<Analysis Complete>>>")
