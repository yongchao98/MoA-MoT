def calculate_glycan_masses():
    """
    Calculates the expected m/z for three modified sialylated glycans.
    The script follows these steps:
    1.  Define monoisotopic masses of elements and chemical groups.
    2.  Calculate the mass of the starting biantennary glycan (A2G2S2).
    3.  Model the chemical reactions (lactonization, amidation) and their effect on mass.
    4.  Model the permethylation and its effect on mass, considering the change in available sites.
    5.  Calculate the final m/z for the singly sodiated ions [M+Na]+.
    """
    # Monoisotopic masses
    mass_H = 1.007825
    mass_C = 12.000000
    mass_O = 15.994915
    mass_N = 14.003074
    mass_Na = 22.989770

    # Base glycan: A2G2S2 -> Man(3)GlcNAc(4)Gal(2)Neu5Ac(2)
    # Molecular Formula: C84 H138 N6 O62
    base_glycan_mass = (84 * mass_C + 138 * mass_H + 6 * mass_N + 62 * mass_O)

    # Mass changes from reactions
    # Lactonization: loss of H2O
    lactone_change = -(2 * mass_H + mass_O)
    # Amidation: -COOH -> -CONH2 (loss of O, gain of NH)
    amide_change = mass_N + mass_H - mass_O
    # Permethylation: addition of CH3 for each H on OH/NH groups, i.e., addition of CH2
    methylation_per_site = mass_C + 2 * mass_H

    # Number of methylation sites on the initial glycan (28 -OH, 6 -NH, 2 -COOH)
    initial_methylation_sites = 36

    glycans = {
        "A2G(4)2S(3)2": {"lactones": 2, "amides": 0},
        "A2G(4)S(3)S(6)": {"lactones": 1, "amides": 1},
        "A2G(4)2S(6)2": {"lactones": 0, "amides": 2},
    }

    print(f"Calculations based on starting glycan A2G2S2 (Man3GlcNAc4Gal2Neu5Ac2), neutral mass: {base_glycan_mass:.4f} Da\n")

    for name, reactions in glycans.items():
        num_lactones = reactions["lactones"]
        num_amides = reactions["amides"]

        # Mass after DMT-MM reaction
        mass_after_reaction = base_glycan_mass + num_lactones * lactone_change + num_amides * amide_change

        # Number of sites for permethylation
        # Each lactone consumes 1 -COOH and 1 -OH -> 2 sites lost
        # Each amide consumes 1 -COOH -> 1 site lost
        methylation_sites = initial_methylation_sites - num_lactones * 2 - num_amides * 1
        
        # Mass added by permethylation
        methylation_mass_add = methylation_sites * methylation_per_site

        # Final mass of the neutral permethylated glycan
        final_neutral_mass = mass_after_reaction + methylation_mass_add
        
        # m/z of the singly sodiated ion [M+Na]+
        final_mz = final_neutral_mass + mass_Na
        
        print(f"For {name}:")
        print(f"This glycan forms {num_lactones} lactone(s) and {num_amides} amide(s).")
        print("The mass calculation is as follows:")
        print(f"  Mass after reaction = Base Mass + Lactone(s) Change + Amide(s) Change")
        print(f"  Mass after reaction = {base_glycan_mass:.4f} + {num_lactones * lactone_change:.4f} + {num_amides * amide_change:.4f} = {mass_after_reaction:.4f} Da")
        print(f"  Number of methylation sites = {initial_methylation_sites} - {num_lactones * 2}(lactone) - {num_amides * 1}(amide) = {methylation_sites}")
        print(f"  Mass from methylation = {methylation_sites} sites * {methylation_per_site:.4f} Da/site = {methylation_mass_add:.4f} Da")
        print(f"  Final [M+Na]+ ion m/z = Mass after reaction + Mass from methylation + Mass of Sodium")
        print(f"  Final m/z = {mass_after_reaction:.4f} + {methylation_mass_add:.4f} + {mass_Na:.4f} = {final_mz:.4f}\n")

# Run the calculation
calculate_glycan_masses()