import pandas as pd
from io import StringIO

def find_calcium_channel_hotspots():
    """
    This function identifies and prints key amino acid residues on the human
    voltage-gated calcium channel beta-1 subunit (CACNB1) that are critical for
    interaction with and modulation of the alpha-1 subunit (Cav2.3).
    The residue positions are based on the canonical sequence of human CACNB1 (UniProt: P54289).
    """

    # Data for residues critical for the physical interaction (binding)
    binding_data = """Residue Position,Amino Acid,Code,Description
186,Leucine,L,"Forms part of the hydrophobic Alpha Interaction Domain (AID) binding pocket."
223,Valine,V,"Key residue in the hydrophobic groove that binds the AID helix."
225,Tyrosine,Y,"Contributes significantly to the binding energy within the hydrophobic pocket."
239,Methionine,M,"Lines the AID binding pocket, making critical hydrophobic contacts."
264,Valine,V,"Essential residue for maintaining the structure of the binding groove."
266,Leucine,L,"Another key hydrophobic residue that directly interacts with the AID."
"""

    # Data for residues critical for fine-tuning gating properties (modulation)
    gating_data = """Residue Position,Amino Acid,Code,Description
140,Tryptophan,W,"Located in the SH3 domain; its mutation alters voltage-dependent inactivation."
163,Aspartate,D,"Located in the Hook region; crucial for modulating channel activation kinetics."
186,Leucine,L,"A dual-role residue; in addition to binding, it's critical for transmitting the modulatory effect on channel activation."
356,Aspartate,D,"Located in the GK domain; involved in allosteric communication that affects channel inactivation."
"""

    binding_df = pd.read_csv(StringIO(binding_data))
    gating_df = pd.read_csv(StringIO(gating_data))

    print("--- 1) Hotspots for Interaction with Alpha-1 Subunit (Binding) ---")
    print("These residues in the beta-1 subunit form the core hydrophobic groove that binds the Alpha Interaction Domain (AID) of the alpha-1 subunit.")
    for index, row in binding_df.iterrows():
        print(f"  - Position: {row['Residue Position']}, Residue: {row['Amino Acid']} ({row['Code']}) - {row['Description']}")

    print("\n--- 2) Hotspots for Gating Modulation of Alpha-1 Subunit ---")
    print("These residues in the beta-1 subunit are critical for fine-tuning the voltage-dependent gating properties of the channel.")
    for index, row in gating_df.iterrows():
        print(f"  - Position: {row['Residue Position']}, Residue: {row['Amino Acid']} ({row['Code']}) - {row['Description']}")

if __name__ == "__main__":
    find_calcium_channel_hotspots()