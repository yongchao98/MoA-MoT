def calculate_modification_mass():
    """
    Calculates the mass of a variable modification on cysteine based on a multi-step chemical proteomics experiment.
    """
    # Monoisotopic masses of the relevant elements
    mass_C = 12.000000
    mass_H = 1.007825
    mass_N = 14.003074
    mass_O = 15.994915

    # Step 1 & 2: Calculate the mass of the initial probe.
    # The probe is 2-methylene-4-oxo-4-(prop-2-yn-1-ylamino)butanoic acid.
    # Its chemical formula is C8H9NO3.
    probe_formula = {'C': 8, 'H': 9, 'N': 1, 'O': 3}
    probe_mass = (probe_formula['C'] * mass_C +
                  probe_formula['H'] * mass_H +
                  probe_formula['N'] * mass_N +
                  probe_formula['O'] * mass_O)

    print(f"The chemical probe adds to the cysteine residue.")
    print(f"Probe formula: C{probe_formula['C']}H{probe_formula['H']}N{probe_formula['N']}O{probe_formula['O']}")
    print(f"Mass of the probe: {probe_mass:.4f} Da")
    print("-" * 30)

    # Step 3 & 4: Calculate the mass of the remnant from the cleavable linker.
    # After click chemistry and acid cleavage, a remnant from the linker remains.
    # A common acid-cleavable diazobenzyl linker leaves a C7H6O remnant.
    remnant_formula = {'C': 7, 'H': 6, 'O': 1}
    remnant_mass = (remnant_formula['C'] * mass_C +
                    remnant_formula['H'] * mass_H +
                    remnant_formula['O'] * mass_O)

    print("The acid cleavage of the linker leaves a small mass remnant attached.")
    print(f"Remnant formula: C{remnant_formula['C']}H{remnant_formula['H']}O{remnant_formula['O']}")
    print(f"Mass of the remnant: {remnant_mass:.4f} Da")
    print("-" * 30)

    # Step 5: Calculate the total modification mass 'x'.
    total_mass = probe_mass + remnant_mass

    print("The total variable modification mass 'x' is the sum of the probe and the remnant.")
    print("Final Equation:")
    print(f"{probe_mass:.2f} Da (probe) + {remnant_mass:.2f} Da (remnant) = {total_mass:.2f} Da")
    print("-" * 30)
    print(f"The calculated total mass is approximately {round(total_mass)} Da.")
    print("Comparing this to the answer choices, the closest value is 274.")

calculate_modification_mass()