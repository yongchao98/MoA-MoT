def get_element_counts(formula):
    """A simple parser for molecular formulas."""
    import re
    elements = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    counts = {}
    for el, count in elements:
        counts[el] = counts.get(el, 0) + int(count) if count else counts.get(el, 0) + 1
    return counts

def calculate_formula_from_components(components):
    """Calculates the final molecular formula from a list of components to add and subtract."""
    final_counts = {}
    for sign, formula in components:
        counts = get_element_counts(formula)
        for el, count in counts.items():
            if sign == '+':
                final_counts[el] = final_counts.get(el, 0) + count
            elif sign == '-':
                final_counts[el] = final_counts.get(el, 0) - count

    # Format the final formula string
    formula_str = ""
    for el in sorted(final_counts.keys()):
        if final_counts[el] > 0:
            formula_str += el
            if final_counts[el] > 1:
                formula_str += str(final_counts[el])
    return formula_str

def main():
    """
    Identifies the structures of products A, B, and C by describing their
    formation and verifying their molecular formulas.
    """
    # Starting Material (SM) formula
    sm_formula = "C9H14N2O2"
    
    # Primary Cycloadduct (P') from Ylide (SM - CO2 - H2) + Methyl Propiolate (MP)
    # Ylide = C8H12N2
    # MP = C4H4O2
    # P' = C12H16N2O2
    p_prime_formula = calculate_formula_from_components([
        ('+', sm_formula), ('-', 'CO2'), ('-', 'H2'), ('+', 'C4H4O2')
    ])

    print("--- Analysis of Reaction Products ---")
    print("-" * 35)

    # --- Product C ---
    # Formed by tautomerization and acetylation of SM.
    # Net reaction: SM + Ac2O - H2O  or SM + C2H2O
    formula_c_calc = calculate_formula_from_components([
        ('+', sm_formula), ('+', 'C2H2O')
    ])
    print("Product C:")
    print("  Proposed Structure: 1-(1-acetyl-2,5-dihydro-1H-pyrrol-2-yl)pyrrolidine-2-carboxylic acid")
    print(f"  Formation: Tautomerization and acetylation of the starting material.")
    print(f"  Given Formula:    C11H16N2O3")
    print(f"  Calculated Formula: {formula_c_calc}")
    print("-" * 35)

    # --- Product B ---
    # Formed from P' by aromatization (-H2) and N-oxidation (+O).
    formula_b_calc = calculate_formula_from_components([
        ('+', p_prime_formula), ('-', 'H2'), ('+', 'O')
    ])
    print("Product B:")
    print("  Proposed Structure: methyl 5-(1-oxy-3,4-dihydro-2H-pyrrol-5-yl)-3H-pyrrolizine-7-carboxylate")
    print(f"  Formation: Aromatization and N-oxidation of the primary cycloadduct.")
    print(f"  Given Formula:    C12H14N2O3")
    print(f"  Calculated Formula: {formula_b_calc}")
    print("-" * 35)

    # --- Product A ---
    # Formed from P' by C-acetylation (+C2H2O) and subsequent reduction of imine (+H2).
    formula_a_calc = calculate_formula_from_components([
        ('+', p_prime_formula), ('+', 'C2H2O'), ('+', 'H2')
    ])
    print("Product A:")
    print("  Proposed Structure: Reduced form of methyl 7-acetyl-5-(3,4-dihydro-2H-pyrrol-5-yl)-2,3,6,7-tetrahydro-1H-pyrrolizine-6-carboxylate")
    print(f"  Formation: C-acetylation of the primary cycloadduct, followed by reduction.")
    print(f"  Given Formula:    C14H20N2O3")
    print(f"  Calculated Formula: {formula_a_calc}")
    print("-" * 35)

if __name__ == "__main__":
    main()