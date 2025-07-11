import sys

def design_liquid_crystal(n):
    """
    Calculates and displays the molecular formula for a single-ring
    liquid crystal based on the specified alkyl chain length 'n'.
    """
    print("Designing a single-ring liquid crystal molecule based on the provided principles.")
    print(f"General Structure: CnH2n+1-Ph-CN (where Ph is a benzene ring)")
    print("-" * 50)

    # Principle F: Start with a specific chain length to target room temperature
    print(f"Using recommended alkyl chain length, n = {n}")

    # Calculate atoms in each component of the molecule
    # Alkyl chain: C_n H_{2n+1}
    c_alkyl = n
    h_alkyl = 2 * n + 1

    # Phenyl ring (para-substituted, C6H4)
    c_phenyl = 6
    h_phenyl = 4

    # Cyano group (CN)
    c_cyano = 1
    n_cyano = 1

    # Sum the atoms for the final molecular formula
    # This fulfills the requirement to "output each number in the final equation"
    total_c = c_alkyl + c_phenyl + c_cyano
    total_h = h_alkyl + h_phenyl
    total_n = n_cyano

    print("\nCalculating the molecular formula step-by-step:")
    print(f"Total Carbon (C) atoms = {c_alkyl} (from alkyl) + {c_phenyl} (from phenyl) + {c_cyano} (from cyano) = {total_c}")
    print(f"Total Hydrogen (H) atoms = {h_alkyl} (from alkyl) + {h_phenyl} (from phenyl) = {total_h}")
    print(f"Total Nitrogen (N) atoms = {n_cyano} (from cyano) = {total_n}")
    print("-" * 50)

    # Display the final resulting molecule and its formula
    molecule_name = f"4-{'pentyl' if n==5 else f'C{n}-alkyl'}benzonitrile"
    print(f"Resulting Molecule: {molecule_name}")
    print(f"Final Molecular Formula: C{total_c}H{total_h}N{total_n}")

if __name__ == '__main__':
    # Based on recommendation F, we start with n=5
    chain_length = 5
    design_liquid_crystal(chain_length)