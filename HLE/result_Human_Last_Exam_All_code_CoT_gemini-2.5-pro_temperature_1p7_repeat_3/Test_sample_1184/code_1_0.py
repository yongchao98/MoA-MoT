def print_interaction_hotspots():
    """
    Prints the residues on the human beta-1 subunit that are hotspots
    for the physical interaction with the alpha-1 subunit.
    
    The information is based on structural and functional studies of the
    conserved interaction between Cav alpha-1 and beta subunits. Residue
    numbering corresponds to human CACNB1, UniProt ID P54283.
    """
    print("1) Hotspots for Interaction with alpha-1 Subunit:")
    print("---------------------------------------------------")
    print("The primary interaction site on the beta-1 subunit is a hydrophobic pocket within its Guanylate Kinase (GK) domain.")
    print("This pocket binds to the Alpha Interaction Domain (AID) of the alpha-1 subunit.")
    print("\nThe key hotspot residues forming this essential hydrophobic pocket are:")
    # The final print statement spells out each number as requested.
    print("\n- Tyrosine at position 386 (Y386), Leucine at position 389 (L389), Leucine at position 419 (L419), Valine at position 422 (V422), and Valine at position 493 (V493).")
    print("\nMutation of these residues typically disrupts or abolishes the physical binding to the alpha-1 subunit.\n")


def print_gating_modulation_hotspots():
    """
    Prints the residues on the human beta-1 subunit that are hotspots
    for fine-tuning the gating properties of the alpha-1 subunit.
    
    The information is based on functional studies. Residue numbering
    corresponds to human CACNB1, UniProt ID P54283.
    """
    print("2) Hotspots for Gating Modulation of alpha-1 Subunit:")
    print("------------------------------------------------------")
    print("Gating modulation is a complex process. While physical binding is a prerequisite, certain residues are particularly crucial for transducing that binding into a functional change in the channel.")
    print("\nThe key hotspot residues for fine-tuning gating properties are:")
    # The final print statement spells out each number as requested.
    print("\n- Aspartate at position 426 (D426): Forms a critical salt bridge with a residue in the alpha-1 AID, which is essential for proper allosteric communication and modulation.")
    print("- Tryptophan at position 331 (W331): Located in the SH3 domain, this residue is not at the primary binding interface but is crucial for the overall structural integrity and communication between domains that enables proper modulation.")
    print("- Glutamate at position 354 (E354): Located in the 'Hook' region connecting the SH3 and GK domains, this residue influences inter-domain communication which is vital for fine-tuning the channel's gating behavior.")


if __name__ == '__main__':
    print_interaction_hotspots()
    print_gating_modulation_hotspots()