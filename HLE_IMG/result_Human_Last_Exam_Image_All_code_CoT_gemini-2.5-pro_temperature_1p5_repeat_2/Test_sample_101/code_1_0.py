import collections

def get_formula_string(atom_counts):
    """Converts a dictionary of atom counts to a standard chemical formula string."""
    formula = ""
    # Standard order: C, then H, then alphabetical for the rest.
    if 'C' in atom_counts and atom_counts['C'] > 0:
        formula += "C" + (str(atom_counts['C']) if atom_counts['C'] > 1 else "")
    if 'H' in atom_counts and atom_counts['H'] > 0:
        formula += "H" + (str(atom_counts['H']) if atom_counts['H'] > 1 else "")
    
    other_atoms = sorted([atom for atom in atom_counts if atom not in ['C', 'H']])
    for atom in other_atoms:
        if atom_counts[atom] > 0:
            formula += atom + (str(atom_counts[atom]) if atom_counts[atom] > 1 else "")
            
    return formula

# Atomic composition of initial reactants
# 3-hydroxy-pyridine-2-carbaldehyde: C6H5NO2
aldehyde = collections.Counter({'C': 6, 'H': 5, 'N': 1, 'O': 2})
# Aniline: C6H7N
aniline = collections.Counter({'C': 6, 'H': 7, 'N': 1})
# Water (lost in condensation): H2O
water = collections.Counter({'H': 2, 'O': 1})
# Hydrogen Cyanide (added in Strecker reaction): HCN
hcn = collections.Counter({'H': 1, 'C': 1, 'N': 1})

# Step 1: Calculate the formula of the imine intermediate.
# Reaction: Aldehyde + Aniline -> Imine + Water
reactants_sum = aldehyde + aniline
imine_intermediate = reactants_sum - water

# Step 2: Calculate the formula of Compound A.
# Reaction: Imine + HCN -> Compound A
compound_A = imine_intermediate + hcn

print("--- Analysis of the Reaction ---")
print("The reaction produces Compound A via a two-step process: imine formation followed by a Strecker-type reaction.")
print("\n--- Identification of Compound A ---")
print("Compound A is an α-aminonitrile with the following properties:")
print("IUPAC Name: 2-(cyano(phenylamino)methyl)pyridin-3-ol")

formula_A = get_formula_string(compound_A)
print(f"Molecular Formula: {formula_A}")

print("\nThe chemical equation for the final step is:")
print(f"Imine ({get_formula_string(imine_intermediate)}) + HCN -> Compound A ({get_formula_string(compound_A)})")

# Explicitly printing the numbers in the final equation as requested
print("\nAtom count for the final equation:")
print(f"Imine (C:{imine_intermediate['C']}, H:{imine_intermediate['H']}, N:{imine_intermediate['N']}, O:{imine_intermediate['O']}) + HCN -> Compound A (C:{compound_A['C']}, H:{compound_A['H']}, N:{compound_A['N']}, O:{compound_A['O']})")
