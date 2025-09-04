import re

class ChemistryAnswerChecker:
    """
    A class to check the logical and factual correctness of the LLM's answer
    to the chemistry problem.
    """

    def __init__(self):
        self.errors = []
        # A simple knowledge base of chemical fragments and their atomic compositions.
        # The composition represents the atoms in the fragment itself.
        self.fragments = {
            # Groups that form the options
            "styryl": {'C': 8, 'H': 7},  # C6H5-CH=CH-
            "4-methylstyryl": {'C': 9, 'H': 9},  # CH3-C6H4-CH=CH-
            "1-phenylprop-1-en-2-yl": {'C': 9, 'H': 9}, # C6H5-CH=C(CH3)-
            "oxirane_terminal": {'C': 2, 'H': 3, 'O': 1}, # -CH(O)CH2
            "oxirane_internal_methyl": {'C': 3, 'H': 5, 'O': 1}, # -CH(O)CH(CH3)-

            # Groups that form the product
            "p-tolyl": {'C': 7, 'H': 7},  # CH3-C6H4-
            "enone_ketone_part": {'C': 4, 'H': 5, 'O': 1},  # -CH=CH-C(=O)-CH3
        }

    def _add_compositions(self, comp1, comp2):
        """Adds the atomic counts of two compositions."""
        result = comp1.copy()
        for atom, count in comp2.items():
            result[atom] = result.get(atom, 0) + count
        return result

    def _get_composition(self, name):
        """Retrieves a fragment's composition from the knowledge base."""
        return self.fragments.get(name, {})

    def _format_formula(self, comp):
        """Formats a composition dictionary into a string like C11H12O1."""
        # Sort by a standard order (C, H, then alphabetical)
        formula = ""
        if 'C' in comp:
            formula += f"C{comp['C']}"
        if 'H' in comp:
            formula += f"H{comp['H']}"
        
        other_atoms = sorted([key for key in comp if key not in ['C', 'H']])
        for atom in other_atoms:
            formula += f"{atom}{comp[atom]}"
        return formula

    def _calculate_dou(self, comp):
        """Calculates the Degree of Unsaturation."""
        # Formula: DoU = C + 1 - H/2 - X/2 + N/2
        c = comp.get('C', 0)
        h = comp.get('H', 0)
        # No halogens (X) or nitrogen (N) in this problem
        return c + 1 - (h / 2)

    def check_answer(self):
        """
        Runs a series of checks to validate the LLM's reasoning.
        """
        # --- Given Information ---
        initial_formula_str = "C11H12O"
        initial_comp = {'C': 11, 'H': 12, 'O': 1}
        
        # --- Check 1: LLM's analysis of the product ---
        llm_product_comp = self._add_compositions(
            self._get_composition("p-tolyl"),
            self._get_composition("enone_ketone_part")
        )
        
        if self._format_formula(llm_product_comp) != self._format_formula(initial_comp):
            self.errors.append(
                f"Constraint Failure: The LLM's proposed product, (E)-4-(p-tolyl)but-3-en-2-one, has a formula of "
                f"{self._format_formula(llm_product_comp)}, which does not match the starting material's formula of {initial_formula_str}."
            )
        
        expected_dou = self._calculate_dou(initial_comp)
        if expected_dou != 6.0:
             self.errors.append(f"Logic Error: The calculated DoU for {initial_formula_str} should be 6, but was {expected_dou}.")
        
        # --- Check 2: LLM's evaluation of the options ---
        
        # Option A: 2-methyl-3-styryloxirane
        comp_A = self._add_compositions(self._get_composition("styryl"), self._get_composition("oxirane_internal_methyl"))
        if "p-tolyl" in self.fragments: # Check if it contains the required group
            if 'phenyl' in "2-methyl-3-styryloxirane": # Simplified check
                pass # Correctly identified as having a phenyl group
            else:
                self.errors.append("Logic Error: LLM correctly rejects A, but the reason (phenyl vs p-tolyl) could be checked more robustly.")

        # Option C: 2-styrylepoxide
        comp_C = self._add_compositions(self._get_composition("styryl"), self._get_composition("oxirane_terminal"))
        if self._format_formula(comp_C) == initial_formula_str:
            self.errors.append(
                f"Constraint Failure: The LLM claims Option C has the wrong formula, but the calculated formula "
                f"({self._format_formula(comp_C)}) matches the starting material."
            )
        elif self._format_formula(comp_C) != "C10H10O1":
             self.errors.append(
                f"Logic Error: The LLM correctly identifies Option C's formula as C10H10O, but my calculation gives {self._format_formula(comp_C)}."
            )

        # Option D: 2-(4-methylstyryl)oxirane (The proposed correct answer)
        comp_D = self._add_compositions(self._get_composition("4-methylstyryl"), self._get_composition("oxirane_terminal"))
        if self._format_formula(comp_D) != initial_formula_str:
            self.errors.append(
                f"Constraint Failure: The final answer, Option D, has a calculated formula of {self._format_formula(comp_D)}, "
                f"which does not match the required formula {initial_formula_str}."
            )
        
        # Check if Option D contains the necessary p-tolyl backbone
        if "4-methylstyryl" not in self.fragments:
             self.errors.append("Logic Error: Cannot verify the structure of Option D.")
        
        # --- Final Verdict ---
        if self.errors:
            print("Incorrect")
            for error in self.errors:
                print(f"- {error}")
        else:
            print("Correct")

# Instantiate the checker and run the validation
checker = ChemistryAnswerChecker()
checker.check_answer()