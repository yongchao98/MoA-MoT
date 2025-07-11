# Protein masses in kDa
protein_masses = {
    'GIKS3': 60,
    'CaPK1': 45,
    'CaPK2': 102,
    'CaPK3': 65,
    'CaPK4': 55,
    'CaPK5': 39,
}

# SEC-MALS observed complex masses in kDa. Use 'None' if no complex peak was reported.
observed_complexes = {
    'GIKS3+CaPK1': None,
    'GIKS3+CaPK2': 222,
    'GIKS3+CaPK3': 185,
    'GIKS3+CaPK4': None,
    'GIKS3+CaPK5': 159,
}

# The control experiment showed GIKS3 runs at 120 kDa, indicating it is a dimer.
giks3_dimer_mass = protein_masses['GIKS3'] * 2
print(f"Analysis of Experiment 1: SEC-MALS Interaction\n")
print(f"GIKS3 monomer mass is {protein_masses['GIKS3']} kDa.")
print(f"The control experiment shows a peak at 120 kDa, so GIKS3 is a dimer in solution.")
print(f"Effective mass of GIKS3 for complex formation is {giks3_dimer_mass} kDa.\n")

print("Calculating expected complex masses and comparing to observations:")
print("-" * 60)

for i in range(1, 6):
    kinase_name = f'CaPK{i}'
    kinase_mass = protein_masses[kinase_name]
    
    # Calculate the expected mass of the complex (GIKS3 dimer + 1 kinase molecule)
    expected_mass = giks3_dimer_mass + kinase_mass
    
    experiment_key = f'GIKS3+{kinase_name}'
    observed_mass = observed_complexes[experiment_key]
    
    print(f"Experiment: GIKS3 + {kinase_name}")
    print(f"GIKS3 dimer mass + {kinase_name} mass = Expected complex mass")
    print(f"{giks3_dimer_mass} kDa + {kinase_mass} kDa = {expected_mass} kDa")
    
    if observed_mass is not None:
        print(f"Observed complex peak: {observed_mass} kDa.")
        if observed_mass == expected_mass:
            print("Result: Interaction confirmed.")
        else:
            print("Result: Discrepancy noted.")
    else:
        print(f"Observed complex peak: Not detected.")
        print("Result: No stable interaction detected.")
    print("-" * 60)

print("\nFinal Conclusion based on all three experiments:")
print("CaPK2, CaPK3, and CaPK4 phosphorylate GIKS3 on Serine 25, leading to its activation.")
print("The SEC-MALS experiment detected stable interaction for CaPK2 and CaPK3, but not for CaPK1 and CaPK4.")
print("Therefore, the statement 'CaPK2 and CaPK3 can phosphorylate GIKS3 on serine 25. CaPK4 does not interact with GIKS3 and CaPK1 does not interact with GIKS3' is a correct, albeit incomplete, summary of the data.")
<<<A>>>