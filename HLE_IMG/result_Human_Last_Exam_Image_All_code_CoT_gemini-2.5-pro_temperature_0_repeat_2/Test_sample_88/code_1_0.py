import re
from collections import Counter

def parse_formula(formula):
    """Parses a chemical formula string into a dictionary of atom counts."""
    atoms = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    atom_counts = Counter()
    for atom, count in atoms:
        atom_counts[atom] += int(count) if count else 1
    return atom_counts

def format_formula(atom_dict):
    """Formats an atom count dictionary back into a formula string."""
    return "".join([f"{atom}{count}" for atom, count in sorted(atom_dict.items())])

def format_counts(atom_dict):
    """Formats an atom count dictionary for printing."""
    return ", ".join([f"{atom}={count}" for atom, count in sorted(atom_dict.items())])

def add_formulas(*dicts):
    """Adds multiple atom count dictionaries."""
    result = Counter()
    for d in dicts:
        result.update(d)
    return result

def subtract_formulas(d1, d2):
    """Subtracts one atom count dictionary from another."""
    result = Counter(d1)
    result.subtract(Counter(d2))
    return result

# --- Define Molecules ---
# Reactants
SM = parse_formula("C9H14N2O2")   # Starting Material: (3,4-dihydro-2H-pyrrol-5-yl)proline
MP = parse_formula("C4H4O2")     # Methyl Propiolate
Ac2O = parse_formula("C4H6O3")    # Acetic Anhydride

# Products
A = parse_formula("C14H20N2O3")
B = parse_formula("C12H14N2O3")
C = parse_formula("C11H16N2O3")

# Byproducts/Side Reagents
AcOH = parse_formula("C2H4O2")   # Acetic Acid
MeOH = parse_formula("CH4O")     # Methanol
CO2 = parse_formula("CO2")       # Carbon Dioxide

# --- Verify Transformations ---

print("Verifying the formation of products A, B, and C.\n")

# 1. Formation of Product C
print("--- Equation for Product C ---")
lhs_c = add_formulas(SM, Ac2O)
rhs_c = add_formulas(C, AcOH)

print(f"Proposed Equation: {format_formula(SM)} + {format_formula(Ac2O)} -> {format_formula(C)} + {format_formula(AcOH)}")
print(f"Reactants side atom counts: {format_counts(lhs_c)}")
print(f"Products side atom counts:  {format_counts(rhs_c)}")
if lhs_c == rhs_c:
    print("Conclusion: The equation is balanced. Product C is the mixed anhydride.\n")
else:
    print("Conclusion: The equation is NOT balanced.\n")


# 2. Formation of Product B
print("--- Equation for Product B ---")
lhs_b = add_formulas(SM, MP)
rhs_b = add_formulas(B, MeOH)

print(f"Proposed Equation: {format_formula(SM)} + {format_formula(MP)} -> {format_formula(B)} + {format_formula(MeOH)}")
print(f"Reactants side atom counts: {format_counts(lhs_b)}")
print(f"Products side atom counts:  {format_counts(rhs_b)}")
if lhs_b == rhs_b:
    print("Conclusion: The equation is balanced.\n")
else:
    print("Conclusion: The equation is NOT balanced.\n")


# 3. Formation of Product A
print("--- Equation for Product A ---")
lhs_a = add_formulas(C, MP)
rhs_a = add_formulas(A, CO2)

print(f"Proposed Equation: {format_formula(C)} + {format_formula(MP)} -> {format_formula(A)} + {format_formula(CO2)}")
print(f"Reactants side atom counts: {format_counts(lhs_a)}")
print(f"Products side atom counts:  {format_counts(rhs_a)}")
if lhs_a == rhs_a:
    print("Conclusion: The equation is balanced.\n")
else:
    print("Conclusion: The equation is NOT balanced.\n")
