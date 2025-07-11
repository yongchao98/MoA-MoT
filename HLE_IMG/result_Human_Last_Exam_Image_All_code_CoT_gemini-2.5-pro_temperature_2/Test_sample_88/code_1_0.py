def format_formula(atoms):
    """Formats a dictionary of atoms into a chemical formula string."""
    formula = ""
    for element in ['C', 'H', 'N', 'O']:
        count = atoms.get(element, 0)
        if count > 0:
            formula += element
            if count > 1:
                formula += str(count)
    return formula

def calculate_formula(operations):
    """
    Performs addition and subtraction on molecular formulas represented as dictionaries.
    operations is a list of tuples, e.g., [('add', formula_dict), ('sub', formula_dict)]
    """
    result = {'C': 0, 'H': 0, 'N': 0, 'O': 0}
    for op, atoms in operations:
        for element, count in atoms.items():
            if op == 'add':
                result[element] += count
            elif op == 'sub':
                result[element] -= count
    return result

# --- Define Molecular Formulas ---
SM = {'C': 9, 'H': 14, 'N': 2, 'O': 2}   # (3,4-dihydro-2H-pyrrol-5-yl)proline
MP = {'C': 4, 'H': 4, 'N': 0, 'O': 2}   # Methyl propiolate
A = {'C': 14, 'H': 20, 'N': 2, 'O': 3}  # Product A
B = {'C': 12, 'H': 14, 'N': 2, 'O': 3}  # Product B
C = {'C': 11, 'H': 16, 'N': 2, 'O': 3}  # Product C

# --- Define other relevant molecules/groups ---
ACETYLATION = {'C': 2, 'H': 2, 'O': 1} # Net addition from acetylation (C2H3O - H)
METHANOL = {'C': 1, 'H': 4, 'O': 1}    # CH3OH
CO2 = {'C': 1, 'O': 2}                 # Carbon Dioxide

# --- Verification for Product C ---
print("--- Analysis for Product C ---")
print(f"Starting Material (SM) formula: {format_formula(SM)}")
print(f"Net change from acetylation (adds C2H3O, removes H): {format_formula(ACETYLATION)}")
result_C_calc = calculate_formula([('add', SM), ('add', ACETYLATION)])
print(f"{format_formula(SM)} + {format_formula(ACETYLATION)} = {format_formula(result_C_calc)}")
print(f"Product C formula: {format_formula(C)}")
print(f"Match: {result_C_calc == C}\n")

# --- Verification for Product B ---
print("--- Analysis for Product B ---")
print(f"Starting Material (SM) formula: {format_formula(SM)}")
print(f"Methyl Propiolate (MP) formula: {format_formula(MP)}")
sm_plus_mp = calculate_formula([('add', SM), ('add', MP)])
print(f"{format_formula(SM)} + {format_formula(MP)} = {format_formula(sm_plus_mp)}")
print(f"Eliminating Methanol (CH3OH): -{format_formula(METHANOL)}")
result_B_calc = calculate_formula([('add', sm_plus_mp), ('sub', METHANOL)])
print(f"Resulting formula: {format_formula(result_B_calc)}")
print(f"Product B formula: {format_formula(B)}")
print(f"Match: {result_B_calc == B}\n")

# --- Verification for Product A ---
print("--- Analysis for Product A ---")
print("This reaction appears to proceed via intermediate C.")
print(f"Intermediate C formula: {format_formula(C)}")
print(f"Methyl Propiolate (MP) formula: {format_formula(MP)}")
c_plus_mp = calculate_formula([('add', C), ('add', MP)])
print(f"{format_formula(C)} + {format_formula(MP)} = {format_formula(c_plus_mp)}")
print(f"Eliminating Carbon Dioxide (CO2): -{format_formula(CO2)}")
result_A_calc = calculate_formula([('add', c_plus_mp), ('sub', CO2)])
print(f"Resulting formula: {format_formula(result_A_calc)}")
print(f"Product A formula: {format_formula(A)}")
print(f"Match: {result_A_calc == A}\n")

# Summary of the formation pathways
print("--- Summary of Deduced Reactions ---")
print(f"Product C: SM ({format_formula(SM)}) + Acetyl_Group-H ({format_formula(ACETYLATION)}) -> C ({format_formula(C)})")
print(f"Product B: SM ({format_formula(SM)}) + MP ({format_formula(MP)}) - MeOH ({format_formula(METHANOL)}) -> B ({format_formula(B)})")
print(f"Product A: C ({format_formula(C)}) + MP ({format_formula(MP)}) - CO2 ({format_formula(CO2)}) -> A ({format_formula(A)})")
