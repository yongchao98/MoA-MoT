from collections import Counter

def format_formula(formula_dict):
    """Formats a formula dictionary into a readable string."""
    return " ".join([f"{elem}={count}" for elem, count in sorted(formula_dict.items())])

def add_formulas(f1, f2):
    """Adds two molecular formulas (represented as Counters)."""
    return f1 + f2

def subtract_formulas(f1, f2):
    """Subtracts one molecular formula from another."""
    return f1 - f2

def print_verification(equation_str, result_formula, target_formula):
    """Prints the verification result for one equation."""
    print(f"Verifying equation: {equation_str}")
    print(f"  - Calculated result formula: {format_formula(result_formula)}")
    print(f"  - Target product formula:    {format_formula(target_formula)}")
    if result_formula == target_formula:
        print("  - Result: The formulas MATCH.\n")
    else:
        print("  - Result: The formulas DO NOT MATCH.\n")

# --- Define Molecular Formulas ---
# Reactants and building blocks
SM = Counter({'C': 9, 'H': 14, 'N': 2, 'O': 2})   # Starting Material
MP = Counter({'C': 4, 'H': 4, 'O': 2})        # Methyl Propiolate
ACETYL = Counter({'C': 2, 'H': 2, 'O': 1})     # Net acetyl group addition (C2H2O)
CO2 = Counter({'C': 1, 'O': 2})              # Carbon Dioxide
CH4O = Counter({'C': 1, 'H': 4, 'O': 1})       # Methanol

# Products
A = Counter({'C': 14, 'H': 20, 'N': 2, 'O': 3})
B = Counter({'C': 12, 'H': 14, 'N': 2, 'O': 3})
C = Counter({'C': 11, 'H': 16, 'N': 2, 'O': 3})

print("--- Verifying the formation of products A, B, and C ---\n")

# --- Verification for Product C ---
# Hypothesis: C is the N-acetylated starting material.
# Equation: C = SM + Acetyl group
c_calc = add_formulas(SM, ACETYL)
print_verification("Product C = Starting Material + Acetyl Group", c_calc, C)
print("The numbers in the equation for C are:")
print(f"C({C['C']})H({C['H']})N({C['N']})O({C['O']}) = C({SM['C']})H({SM['H']})N({SM['N']})O({SM['O']}) + C({ACETYL['C']})H({ACETYL['H']})O({ACETYL['O']})\n")


# --- Verification for Product A ---
# Hypothesis: A is formed from C and methyl propiolate, with loss of CO2.
# Equation: A = C + MP - CO2
a_calc = subtract_formulas(add_formulas(C, MP), CO2)
print_verification("Product A = Product C + Methyl Propiolate - CO2", a_calc, A)
print("The numbers in the equation for A are:")
print(f"C({A['C']})H({A['H']})N({A['N']})O({A['O']}) = C({C['C']})H({C['H']})N({C['N']})O({C['O']}) + C({MP['C']})H({MP['H']})O({MP['O']}) - C({CO2['C']})O({CO2['O']})\n")

# --- Verification for Product B ---
# Hypothesis: B is formed from the starting material and methyl propiolate, with loss of methanol.
# Equation: B = SM + MP - CH4O
b_calc = subtract_formulas(add_formulas(SM, MP), CH4O)
print_verification("Product B = Starting Material + Methyl Propiolate - Methanol", b_calc, B)
print("The numbers in the equation for B are:")
print(f"C({B['C']})H({B['H']})N({B['N']})O({B['O']}) = C({SM['C']})H({SM['H']})N({SM['N']})O({SM['O']}) + C({MP['C']})H({MP['H']})O({MP['O']}) - C({CH4O['C']})H({CH4O['H']})O({CH4O['O']})\n")
