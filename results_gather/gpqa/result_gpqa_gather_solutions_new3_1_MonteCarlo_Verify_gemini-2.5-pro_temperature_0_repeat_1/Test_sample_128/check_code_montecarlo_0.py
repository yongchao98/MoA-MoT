import re

class ChemistryProblemChecker:
    """
    A class to verify the solution to a multi-step organic chemistry problem.
    It simulates the reaction pathway using string manipulation based on established chemical rules.
    """

    def __init__(self, question_data):
        """
        Initializes the checker with the problem's data.
        
        Args:
            question_data (dict): A dictionary containing the options, hints, and the LLM's answer key.
        """
        self.options = question_data["options"]
        self.hints = question_data["hints"]
        self.llm_answer_key = question_data["llm_answer_key"]
        self.errors = []
        self.log = []

    def _add_error(self, message):
        """Adds an error message to the list of errors."""
        self.errors.append(message)

    def _log_step(self, message):
        """Logs a successful step in the reasoning process."""
        self.log.append(message)

    def check_compound_a(self):
        """
        Step 1: Deduce Compound A from the Wittig reaction (Hint a) and verify with IR (Hint b).
        """
        wittig_product = self.hints['a']
        # Perform a retro-Wittig analysis on the product: 1,2-dimethyl-4-(propan-2-ylidene)cyclopentane.
        # The ylide part is =C(CH3)2. Replacing it with =O gives a ketone at position 4.
        # The ketone is 1,2-dimethylcyclopentan-4-one.
        # According to IUPAC rules, the carbonyl gets C1, so renumbering gives 3,4-dimethylcyclopentan-1-one.
        deduced_a = "3,4-dimethylcyclopentan-1-one"
        self._log_step(f"Step 1 (Deduce A): Based on retro-Wittig analysis of '{wittig_product}', Compound A is determined to be '{deduced_a}'.")

        # Verify with IR hint for A
        ir_a = self.hints['b']['A']
        if "cyclopentan" in deduced_a and "1750" in ir_a:
            self._log_step("Step 1 (Verify A): The deduced structure is a cyclopentanone, which is consistent with the characteristic IR peak at ~1750 cm^-1 due to ring strain.")
            return deduced_a
        else:
            self._add_error(f"Constraint Mismatch (Hint b): Deduced Compound A ('{deduced_a}') is not consistent with its IR hint ('{ir_a}'). A cyclopentanone is expected for a ~1750 cm^-1 peak.")
            return None

    def check_reaction_a_to_c(self, compound_a):
        """
        Step 2 & 3: Trace the reactions from A to B (cyanohydrin) and B to C (nitrile reduction).
        """
        if not compound_a:
            return None
        # A + HCN -> B (cyanohydrin formation)
        # 3,4-dimethylcyclopentan-1-one -> 1-cyano-3,4-dimethylcyclopentan-1-ol
        deduced_b = "1-cyano-3,4-dimethylcyclopentan-1-ol"
        self._log_step(f"Step 2 (A -> B): Reaction of Compound A with HCN yields the cyanohydrin '{deduced_b}'.")

        # B + H2/Pd -> C (nitrile reduction)
        # 1-cyano... -> 1-(aminomethyl)...
        deduced_c = "1-(aminomethyl)-3,4-dimethylcyclopentan-1-ol"
        self._log_step(f"Step 3 (B -> C): Reduction of the nitrile in Compound B yields the amino alcohol '{deduced_c}'.")
        return deduced_c

    def check_reaction_c_to_e(self, compound_c):
        """
        Step 4: Trace the Tiffeneau-Demjanov rearrangement from C to E.
        """
        if not compound_c:
            return None
        # C + HNO2 -> D -> E (ring expansion)
        # This is a Tiffeneau-Demjanov rearrangement.
        # Key transformations:
        # - Ring size: cyclopentane -> cyclohexane
        # - Functional group: 1-(aminomethyl)...-1-ol -> 1-one
        # - Substituents: 3,4-dimethyl are preserved.
        if "1-(aminomethyl)-3,4-dimethylcyclopentan-1-ol" in compound_c:
            deduced_e = "3,4-dimethylcyclohexan-1-one"
            self._log_step(f"Step 4 (C -> E): The Tiffeneau-Demjanov rearrangement of Compound C results in a ring expansion to yield '{deduced_e}'.")
            
            # Verify with IR hint for E
            ir_e = self.hints['b']['E']
            if "cyclohexan" in deduced_e and "1715" in ir_e:
                self._log_step("Step 4 (Verify E): The deduced structure is a cyclohexanone, which is consistent with the characteristic IR peak at ~1715 cm^-1.")
                return deduced_e
            else:
                self._add_error(f"Constraint Mismatch (Hint b): Deduced Compound E ('{deduced_e}') is not consistent with its IR hint ('{ir_e}'). A cyclohexanone is expected for a ~1715 cm^-1 peak.")
                return None
        else:
            self._add_error("Step 4 Error: The input compound for the rearrangement step is incorrect, breaking the logical chain.")
            return None

    def run_check(self):
        """
        Executes the full verification process and returns the final result.
        """
        compound_a = self.check_compound_a()
        compound_c = self.check_reaction_a_to_c(compound_a)
        deduced_e = self.check_reaction_c_to_e(compound_c)

        if not deduced_e:
            # Errors have already been logged by the failing step.
            return f"Incorrect. Reasons:\n" + "\n".join(self.errors)

        # Final verification against options and LLM's choice
        correct_option_key = None
        for key, value in self.options.items():
            if value == deduced_e:
                correct_option_key = key
                break
        
        if not correct_option_key:
            self._add_error(f"Final Product Error: The logically deduced final product '{deduced_e}' is not present in the given options.")
        else:
            self._log_step(f"Step 5 (Match Option): Deduced product '{deduced_e}' correctly matches option {correct_option_key}.")
            if self.llm_answer_key != correct_option_key:
                self._add_error(f"Answer Choice Mismatch: The correct option is {correct_option_key}, but the LLM chose {self.llm_answer_key}.")

        if not self.errors:
            return "Correct"
        else:
            return f"Incorrect. Reasons:\n" + "\n".join(self.errors)

# Data extracted from the problem description
question_data = {
    "options": {
        "A": "2,2,3,4-tetramethylcyclobutan-1-one",
        "B": "4-methylcycloheptan-1-one",
        "C": "2,3,4-trimethylcyclopentan-1-one",
        "D": "3,4-dimethylcyclohexan-1-one"
    },
    "hints": {
        "a": "1,2-dimethyl-4-(propan-2-ylidene)cyclopentane",
        "b": {
            "A": "strong peak at around 1750 cm^-1",
            "E": "strong peak at around 1715 cm^-1"
        }
    },
    "llm_answer_key": "D" # Extracted from "<<<D>>>"
}

# Instantiate the checker and run the verification
checker = ChemistryProblemChecker(question_data)
result = checker.run_check()
print(result)