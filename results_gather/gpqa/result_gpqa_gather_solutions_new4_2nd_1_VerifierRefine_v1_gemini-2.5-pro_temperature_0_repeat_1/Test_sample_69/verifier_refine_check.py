import sys
from io import StringIO

class ChemicalKnowledgeBase:
    """A simple knowledge base for chemical properties and reactions based on the puzzle's logic."""
    def __init__(self):
        self.properties = {
            'S8': {'state': 'solid', 'name': 'Sulfur (S₈)'},
            'Cl2': {'state': 'gas', 'name': 'Chlorine (Cl₂)'},
            'SCl2': {'state': 'liquid', 'color': 'red', 'name': 'Sulfur dichloride (SCl₂)'},
            'C2H4': {'state': 'gas', 'name': 'Ethene (C₂H₄)'},
            'H2O': {'state': 'liquid', 'name': 'Water (H₂O)'},
            'HCl': {'acid_strength': 'strong', 'name': 'Hydrochloric acid (HCl)'},
            'H2SO3': {'acid_strength': 'weak', 'name': 'Sulfurous acid (H₂SO₃)'},
            'C2H4Cl2': {'use': 'solvent', 'name': '1,2-dichloroethane (C₂H₄Cl₂)'},
            'mustard_gas': {
                'formula': '(ClCH2CH2)2S',
                'hazard': 'extremely hazardous',
                'symmetry': 'C2',
                'name': 'Mustard gas'
            },
        }

    def get_property(self, substance, prop):
        return self.properties.get(substance, {}).get(prop)

    def check_reaction(self, reaction_name, reactants, products, stoichiometry):
        """Checks if a reaction is valid based on the puzzle's most consistent interpretation."""
        # Reaction 1: A(s) + 8 B(g) -> C
        if reaction_name == "reaction_1":
            A, B = reactants
            C, = products
            ratio_A, ratio_B = stoichiometry
            # The reaction S₈ + 8Cl₂ → 8SCl₂ fits the 1:8 ratio.
            return (A == 'S8' and B == 'Cl2' and C == 'SCl2' and
                    ratio_A == 1 and ratio_B == 8)

        # Reaction 3: C + H2O -> A + F + G
        if reaction_name == "reaction_3":
            C, H2O = reactants
            A, F, G = products
            # The hydrolysis of SCl₂ produces S, HCl, and SO₂ (which forms H₂SO₃)
            return (C == 'SCl2' and H2O == 'H2O' and A == 'S8' and
                    F == 'HCl' and G == 'H2SO3')

        # Reaction 4: D(g) + B(g) -> H
        if reaction_name == "reaction_4":
            D, B = reactants
            H, = products
            ratio_D, ratio_B = stoichiometry
            # The reaction C₂H₄ + Cl₂ → C₂H₄Cl₂ fits the 1:1 ratio.
            return (D == 'C2H4' and B == 'Cl2' and H == 'C2H4Cl2' and
                    ratio_D == 1 and ratio_B == 1)

        # Reaction 2: C + 2 D(g) -> E
        if reaction_name == "reaction_2":
            C, D = reactants
            E, = products
            ratio_C, ratio_D = stoichiometry
            # The Levinstein process: SCl₂ + 2C₂H₄ → (ClCH₂CH₂)₂S
            return (C == 'SCl2' and D == 'C2H4' and E == 'mustard_gas' and
                    ratio_C == 1 and ratio_D == 2)

        return False

def check_correctness():
    """
    Checks the correctness of the final answer by verifying each step of the deduction.
    """
    try:
        # 1. Define the proposed solution from the final answer text
        proposed_solution = {
            'A': 'S8',
            'B': 'Cl2',
            'C': 'SCl2',
            'D': 'C2H4',
            'E': 'mustard_gas',
            'F': 'HCl',
            'G': 'H2SO3',
            'H': 'C2H4Cl2',
            'final_symmetry': 'C2',
            'final_option': 'A'
        }
        
        # 2. Initialize the knowledge base
        kb = ChemicalKnowledgeBase()
        
        # 3. Perform checks for each clue
        
        # Check Clue 1: A(s) + 8 B(g) -> C (bright red product)
        A, B, C = proposed_solution['A'], proposed_solution['B'], proposed_solution['C']
        if kb.get_property(A, 'state') != 'solid':
            return f"Constraint check failed for Clue 1: Proposed A ({kb.get_property(A, 'name')}) is not a solid."
        if kb.get_property(B, 'state') != 'gas':
            return f"Constraint check failed for Clue 1: Proposed B ({kb.get_property(B, 'name')}) is not a gas."
        if not kb.check_reaction("reaction_1", [A, B], [C], (1, 8)):
            return f"Constraint check failed for Clue 1: The reaction {A} + 8{B} -> {C} is not consistent with the 1:8 stoichiometry."
        if kb.get_property(C, 'color') != 'red':
            return f"Constraint check failed for Clue 1: Proposed C ({kb.get_property(C, 'name')}) is not a red product."

        # Check Clue 3: C + H2O -> A(s) + F(strong acid) + G(weak acid)
        F, G = proposed_solution['F'], proposed_solution['G']
        if not kb.check_reaction("reaction_3", [C, 'H2O'], [A, F, G], None):
            return f"Constraint check failed for Clue 3: The hydrolysis of {C} does not produce {A}, {F}, and {G} as described."
        if kb.get_property(F, 'acid_strength') != 'strong':
            return f"Constraint check failed for Clue 3: Proposed F ({kb.get_property(F, 'name')}) is not a strong acid."
        if kb.get_property(G, 'acid_strength') != 'weak':
            return f"Constraint check failed for Clue 3: Proposed G ({kb.get_property(G, 'name')}) is not a weak acid."

        # Check Clue 4: D(g) + B(g) -> H (solvent) (1:1 ratio)
        D, H = proposed_solution['D'], proposed_solution['H']
        if kb.get_property(D, 'state') != 'gas':
            return f"Constraint check failed for Clue 4: Proposed D ({kb.get_property(D, 'name')}) is not a gas."
        if not kb.check_reaction("reaction_4", [D, B], [H], (1, 1)):
            return f"Constraint check failed for Clue 4: The reaction {D} + {B} -> {H} is not consistent with the 1:1 stoichiometry."
        if kb.get_property(H, 'use') != 'solvent':
            return f"Constraint check failed for Clue 4: Proposed H ({kb.get_property(H, 'name')}) is not a known solvent."

        # Check Clue 2: C + 2 D(g) -> E (extremely hazardous product)
        E = proposed_solution['E']
        if not kb.check_reaction("reaction_2", [C, D], [E], (1, 2)):
            return f"Constraint check failed for Clue 2: The reaction {C} + 2{D} -> {E} is not consistent with the 1:2 stoichiometry."
        if kb.get_property(E, 'hazard') != 'extremely hazardous':
            return f"Constraint check failed for Clue 2: Proposed E ({kb.get_property(E, 'name')}) does not fit the 'extremely hazardous' description."

        # Check Final Question: What is the molecular symmetry group of E?
        expected_symmetry = kb.get_property(E, 'symmetry')
        if proposed_solution['final_symmetry'] != expected_symmetry:
            return f"Final Answer check failed: The proposed symmetry ({proposed_solution['final_symmetry']}) is incorrect for E ({kb.get_property(E, 'name')}), which has {expected_symmetry} symmetry."

        # Check if the final option matches the conclusion
        # Options: A) C2, B) C2v, C) D∞h, D) D4h
        options = {'A': 'C2', 'B': 'C2v', 'C': 'D∞h', 'D': 'D4h'}
        if options.get(proposed_solution['final_option']) != proposed_solution['final_symmetry']:
            return f"Final Answer check failed: The final option '{proposed_solution['final_option']}' does not correspond to the derived symmetry '{proposed_solution['final_symmetry']}'."

        return "Correct"
    except Exception as e:
        return f"An error occurred during the check: {e}"

# Execute the check and print the result
result = check_correctness()
print(result)