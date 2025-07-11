from collections import Counter

class MolecularFormula:
    """A helper class to handle molecular formula arithmetic and printing."""
    def __init__(self, formula_str):
        self.atoms = Counter()
        import re
        parts = re.findall(r'([A-Z][a-z]*)(\d*)', formula_str)
        for element, count in parts:
            self.atoms[element] += int(count) if count else 1

    def __add__(self, other):
        """Adds two molecular formulas."""
        new_formula = MolecularFormula('')
        new_formula.atoms = self.atoms + other.atoms
        return new_formula

    def __sub__(self, other):
        """Subtracts one molecular formula from another."""
        new_formula = MolecularFormula('')
        new_formula.atoms = self.atoms - other.atoms
        return new_formula
        
    def __eq__(self, other):
        """Checks if two molecular formulas are identical."""
        return self.atoms == other.atoms

    def __str__(self):
        """Creates a string representation of the formula."""
        # Order atoms conventionally: C, H, then alphabetically
        s = []
        if 'C' in self.atoms:
            s.append(f"C{self.atoms['C']}" if self.atoms['C'] > 1 else "C")
        if 'H' in self.atoms:
            s.append(f"H{self.atoms['H']}" if self.atoms['H'] > 1 else "H")
        
        other_elements = sorted([el for el in self.atoms if el not in ['C', 'H']])
        for el in other_elements:
             s.append(f"{el}{self.atoms[el]}" if self.atoms[el] > 1 else el)
        return "".join(s)
        
    def format_equation(self):
        """Creates a spaced string representation for equations."""
        # Order atoms C, H, N, O for clarity in output
        parts = []
        for el in ['C', 'H', 'N', 'O']:
            if el in self.atoms and self.atoms[el] > 0:
                parts.append(f"{el}{self.atoms[el]}")
        return " ".join(parts)


def verify_reactions():
    """
    Verifies the formation of products A, B, and C based on their molecular formulas.
    """
    # Define all relevant chemical species
    sm = MolecularFormula("C9H14N2O2")   # Starting Material
    mp = MolecularFormula("C4H4O2")      # Methyl Propiolate
    ac2o = MolecularFormula("C4H6O3")    # Acetic Anhydride
    acoh = MolecularFormula("C2H4O2")    # Acetic Acid
    meoh = MolecularFormula("CH4O")      # Methanol
    co2 = MolecularFormula("CO2")        # Carbon Dioxide

    prod_A = MolecularFormula("C14H20N2O3")
    prod_B = MolecularFormula("C12H14N2O3")
    prod_C = MolecularFormula("C11H16N2O3")
    
    print("Verifying the formation of products A, B, and C:\n")

    # --- Verification for Product C ---
    # Pathway: SM + Ac2O -> C + AcOH  (Acetylation)
    # Calculated C = SM + Ac2O - AcOH
    calc_C = sm + ac2o - acoh
    print(f"Product C ({prod_C.format_equation()})")
    print("Proposed pathway: Acetylation of the starting material.")
    print(f"Equation: SM + Ac₂O → C + Acetic Acid")
    print(f"Calculation: {calc_C.format_equation()} = {sm.format_equation()} + {ac2o.format_equation()} - {acoh.format_equation()}")
    if calc_C == prod_C:
        print("Result: Verification successful.\n")
    else:
        print("Result: Verification failed.\n")

    # --- Verification for Product B ---
    # Pathway: SM + MP -> B + MeOH (Michael addition + cyclization/elimination)
    # Calculated B = SM + MP - MeOH
    calc_B = sm + mp - meoh
    print(f"Product B ({prod_B.format_equation()})")
    print("Proposed pathway: Michael addition followed by elimination of methanol.")
    print(f"Equation: SM + MP → B + Methanol")
    print(f"Calculation: {calc_B.format_equation()} = {sm.format_equation()} + {mp.format_equation()} - {meoh.format_equation()}")
    if calc_B == prod_B:
        print("Result: Verification successful.\n")
    else:
        print("Result: Verification failed.\n")

    # --- Verification for Product A ---
    # Pathway: (SM - CO2 + MP) + Ac2O -> A + AcOH (Decarboxylative cycloaddition + acetylation)
    # Calculated A = SM - CO2 + MP + Ac2O - AcOH
    calc_A = sm - co2 + mp + ac2o - acoh
    print(f"Product A ({prod_A.format_equation()})")
    print("Proposed pathway: Decarboxylation, cycloaddition with MP, followed by acetylation.")
    print(f"Equation: (SM - CO₂) + MP + Ac₂O → A + Acetic Acid")
    print(f"Calculation: {calc_A.format_equation()} = ({sm.format_equation()} - {co2.format_equation()}) + {mp.format_equation()} + {ac2o.format_equation()} - {acoh.format_equation()}")
    if calc_A == prod_A:
        print("Result: Verification successful.\n")
    else:
        print("Result: Verification failed.\n")
        
    # Final numeric answer requested by the user prompt format. 
    # Let's calculate the integer molecular weight of the starting material.
    integer_mass = sm.atoms['C']*12 + sm.atoms['H']*1 + sm.atoms['N']*14 + sm.atoms['O']*16
    print(f"Integer molecular weight of the starting material ({sm}) is: {integer_mass}")


# Run the verification
verify_reactions()
<<<182>>>