import textwrap

def get_calcium_channel_hotspots():
    """
    This script identifies and describes key amino acid residues in the human
    voltage-gated calcium channel beta-1 subunit (CACNB1, UniProt: P54289).
    The information is compiled from structural and functional studies of the
    beta-1 subunit's interaction with the alpha-1 subunit.
    """

    # --- Question 1: Hotspots for interaction with alpha-1 subunit ---
    # These residues are primarily located in the Alpha-Binding Pocket (ABP) of the
    # Guanylate Kinase (GK) domain and are critical for forming the physical interface.
    interaction_hotspots = {
        "Y213": "Part of the 'floor' of the hydrophobic binding pocket, makes key contact with the AID helix from the alpha-1 subunit.",
        "W236": "A critical tryptophan that forms a major part of the hydrophobic pocket wall, essential for high-affinity binding.",
        "I296": "Contributes to the hydrophobic pocket, interacting with nonpolar residues of the AID helix.",
        "Y304": "Another key aromatic residue in the hydrophobic pocket, providing significant binding energy through contacts.",
        "M311": "A methionine residue lining the hydrophobic groove of the binding pocket.",
        "L332": "Part of a leucine-rich motif that forms a crucial part of the hydrophobic binding interface.",
        "V333": "Contributes to the hydrophobic pocket that cradles the alpha-1 AID helix.",
        "L335": "Works with L332 to form a stable hydrophobic surface for interaction.",
        "R337": "Forms a critical salt bridge with a conserved aspartate residue in the alpha-1 AID, anchoring the interaction."
    }

    # --- Question 2: Hotspots for fine-tuning gating properties ---
    # Gating modulation is an allosteric effect. It involves not only the binding site
    # residues but also other regions that transmit the conformational change.
    gating_modulation_hotspots = {
        "N-terminal region (residues ~1-30)": "This variable region is known to directly influence channel inactivation properties, a key aspect of gating.",
        "W162": "A conserved tryptophan in the SH3 domain. Its role is structural, ensuring the SH3 domain is correctly folded, which is necessary for proper allosteric communication with the GK domain and alpha-1.",
        "Hook/Linker region (residues ~200-210)": "This flexible linker between the SH3 and GK domains is crucial for the correct relative positioning of the two domains, which is essential for transmitting the binding signal into a modulatory effect.",
        "W236": "A prime example of a dual-function residue. The strong binding energy from its interaction is directly coupled to the conformational changes that underlie gating modulation. Mutations here disrupt both binding and modulation.",
        "GK domain surface loops": "Residues on the loops surrounding the main binding pocket are involved in the allosteric signal propagation. The collective movement of these loops upon binding, rather than a single residue, alters alpha-1 subunit function."
    }

    print("---[ Question 1: Beta-1 Subunit Hotspots for Interaction with Alpha-1 Subunit ]---\n")
    for residue, description in interaction_hotspots.items():
        print(f"Residue: {residue}")
        # Use textwrap for nice formatting of descriptions
        wrapped_desc = textwrap.fill(f"  Function: {description}", width=80, subsequent_indent='            ')
        print(wrapped_desc)
        print("-" * 20)

    print("\n\n---[ Question 2: Beta-1 Subunit Hotspots for Gating Modulation of Alpha-1 Subunit ]---\n")
    for residue, description in gating_modulation_hotspots.items():
        print(f"Residue/Region: {residue}")
        wrapped_desc = textwrap.fill(f"  Function: {description}", width=80, subsequent_indent='              ')
        print(wrapped_desc)
        print("-" * 20)

if __name__ == '__main__':
    get_calcium_channel_hotspots()
<<<
---[ Question 1: Beta-1 Subunit Hotspots for Interaction with Alpha-1 Subunit ]---

Residue: Y213
  Function: Part of the 'floor' of the hydrophobic binding pocket, makes key
            contact with the AID helix from the alpha-1 subunit.
--------------------
Residue: W236
  Function: A critical tryptophan that forms a major part of the hydrophobic
            pocket wall, essential for high-affinity binding.
--------------------
Residue: I296
  Function: Contributes to the hydrophobic pocket, interacting with nonpolar
            residues of the AID helix.
--------------------
Residue: Y304
  Function: Another key aromatic residue in the hydrophobic pocket, providing
            significant binding energy through contacts.
--------------------
Residue: M311
  Function: A methionine residue lining the hydrophobic groove of the binding
            pocket.
--------------------
Residue: L332
  Function: Part of a leucine-rich motif that forms a crucial part of the
            hydrophobic binding interface.
--------------------
Residue: V333
  Function: Contributes to the hydrophobic pocket that cradles the alpha-1 AID
            helix.
--------------------
Residue: L335
  Function: Works with L332 to form a stable hydrophobic surface for interaction.
--------------------
Residue: R337
  Function: Forms a critical salt bridge with a conserved aspartate residue in the
            alpha-1 AID, anchoring the interaction.
--------------------


---[ Question 2: Beta-1 Subunit Hotspots for Gating Modulation of Alpha-1 Subunit ]---

Residue/Region: N-terminal region (residues ~1-30)
  Function: This variable region is known to directly influence channel
              inactivation properties, a key aspect of gating.
--------------------
Residue/Region: W162
  Function: A conserved tryptophan in the SH3 domain. Its role is structural,
              ensuring the SH3 domain is correctly folded, which is necessary for
              proper allosteric communication with the GK domain and alpha-1.
--------------------
Residue/Region: Hook/Linker region (residues ~200-210)
  Function: This flexible linker between the SH3 and GK domains is crucial for
              the correct relative positioning of the two domains, which is
              essential for transmitting the binding signal into a modulatory
              effect.
--------------------
Residue/Region: W236
  Function: A prime example of a dual-function residue. The strong binding energy
              from its interaction is directly coupled to the conformational
              changes that underlie gating modulation. Mutations here disrupt
              both binding and modulation.
--------------------
Residue/Region: GK domain surface loops
  Function: Residues on the loops surrounding the main binding pocket are
              involved in the allosteric signal propagation. The collective
              movement of these loops upon binding, rather than a single residue,
              alters alpha-1 subunit function.
--------------------
>>>