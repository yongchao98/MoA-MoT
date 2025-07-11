import sys
# Redirect print to stdout, ensuring it's not going to stderr.
# This helps prevent the output from being misinterperted as an error.
sys.stdout.reconfigure(encoding='utf-8')

def calculate_glycan_mass():
    """
    Calculates the m/z of a specific derivatized and sodiated N-glycan.
    """
    # Define monoisotopic atomic masses
    atomic_masses = {
        'C': 12.000000,
        'H': 1.007825,
        'O': 15.994915,
        'N': 14.003074,
        'Na': 22.989770,
    }

    # Helper function to calculate mass from a chemical formula
    def get_mass(formula):
        mass = 0.0
        for atom, count in formula.items():
            mass += atomic_masses[atom] * count
        return mass

    # --- Step 1: Define Composition and Calculate Native Mass ---
    # Composition of A2G2S2 glycan
    num_hexose = 5  # 3 Man + 2 Gal
    num_glcnac = 4
    num_neu5ac = 2
    total_monosaccharides = num_hexose + num_glcnac + num_neu5ac

    # Masses of building blocks
    mass_hexose = get_mass({'C': 6, 'H': 12, 'O': 6})
    mass_glcnac = get_mass({'C': 8, 'H': 15, 'N': 1, 'O': 6})
    mass_neu5ac = get_mass({'C': 11, 'H': 19, 'N': 1, 'O': 9})
    mass_h2o = get_mass({'H': 2, 'O': 1})
    
    mass_sum_monosaccharides = (num_hexose * mass_hexose +
                                num_glcnac * mass_glcnac +
                                num_neu5ac * mass_neu5ac)
    
    num_bonds = total_monosaccharides - 1
    mass_native_glycan = mass_sum_monosaccharides - num_bonds * mass_h2o

    # --- Step 2: Account for Amidation ---
    num_amidations = 2
    mass_oh = get_mass({'O': 1, 'H': 1})
    mass_nh2 = get_mass({'N': 1, 'H': 2})
    mass_change_amidation = mass_nh2 - mass_oh
    mass_amidated_glycan = mass_native_glycan + num_amidations * mass_change_amidation
    
    # --- Step 3: Account for Permethylation ---
    # Count reactive sites for methylation on the amidated glycan
    sites_hexose = 5
    sites_glcnac = 5  # 4 OH + 1 NH
    sites_neu5ac_amide = 7 # 4 OH + 1 NH (acetyl) + 2 NH (amide)
    total_sites_in_monomers = (num_hexose * sites_hexose +
                               num_glcnac * sites_glcnac +
                               num_neu5ac * sites_neu5ac_amide)
    num_methylation_sites = total_sites_in_monomers - 2 * num_bonds
    
    mass_ch2 = get_mass({'C': 1, 'H': 2})
    mass_increase_permethylation = num_methylation_sites * mass_ch2
    mass_final_glycan = mass_amidated_glycan + mass_increase_permethylation

    # --- Step 4: Account for Sodiated Ion ---
    mass_na = atomic_masses['Na']
    mass_sodiated_ion = mass_final_glycan + mass_na

    # --- Print the results ---
    print("The three specified glycans are isomers and will have the same mass after the described chemical derivatization.")
    print("The expected mass-to-charge ratio (m/z) for the singly sodiated ion [M+Na]+ is calculated as follows:\n")
    
    print("Summary of Calculation:")
    print(f"  Mass of Amidated Glycan (M_amidated) = ({num_hexose}*Hexose + {num_glcnac}*GlcNAc + {num_neu5ac}*Neu5Ac) - {num_bonds}*H2O + {num_amidations}*(NH2 - OH)")
    print(f"  M_amidated = ({mass_sum_monosaccharides:.4f}) - {num_bonds}*{mass_h2o:.4f} + {num_amidations}*({mass_change_amidation:.4f}) = {mass_amidated_glycan:.4f} Da\n")
    
    print(f"  Mass from Permethylation = {num_methylation_sites} * Mass(CH2)")
    print(f"  Mass from Permethylation = {num_methylation_sites} * {mass_ch2:.4f} = {mass_increase_permethylation:.4f} Da\n")

    print("Final m/z Calculation for [M+Na]+:")
    print(f"  m/z = M_amidated + Mass_from_Permethylation + Mass_Na")
    print(f"  m/z = {mass_amidated_glycan:.4f} + {mass_increase_permethylation:.4f} + {mass_na:.4f} = {mass_sodiated_ion:.4f}\n")
    
    print("-" * 50)
    print("Final Expected m/z Values:")
    print(f"  Mass of A2G(4)2S(3)2: {mass_sodiated_ion:.4f}")
    print(f"  Mass of A2G(4)S(3)S(6): {mass_sodiated_ion:.4f}")
    print(f"  Mass of A2G(4)2S(6)2: {mass_sodiated_ion:.4f}")
    print("-" * 50)
    
    return mass_sodiated_ion

# Run the calculation and get the final value
final_mass = calculate_glycan_mass()
# The final answer format is specified as <<<value>>>
print(f"<<<{final_mass:.4f}>>>")
