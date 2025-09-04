import re

class ChemistryChecker:
    """
    A class to verify the logical steps in a multi-step organic chemistry problem.
    It focuses on the most complex and error-prone part: the Diels-Alder stereochemistry.
    """

    def __init__(self, chosen_answer, reasoning):
        self.chosen_answer = chosen_answer
        self.reasoning = reasoning
        self.errors = []
        self.options = {
            "A": "(1S,4S)-4,6,6-trimethylcyclohex-2-enol",
            "B": "(1S,4R,5S,6R)-4,5,6-trimethylcyclohex-2-enol",
            "C": "(1S,4R,5S,6S)-4,5,6-trimethylcyclohex-2-enol",
            "D": "(1S,4R)-4,6,6-trimethylcyclohex-2-enol"
        }

    def get_relative_stereochem_from_rs(self, config1, config2):
        """
        Determines relative stereochemistry (cis/trans) for adjacent carbons on a 
        cyclohexane ring from their R/S labels.
        Rule: (R,S) or (S,R) is cis. (R,R) or (S,S) is trans.
        """
        c1 = config1[0].upper()
        c2 = config2[0].upper()
        return "cis" if c1 != c2 else "trans"

    def check_product_skeleton(self):
        """
        Verifies the basic structure (skeleton) of the Diels-Alder product.
        """
        correct_skeleton = "4,5,6-trimethylcyclohex-2-enol"
        chosen_name = self.options.get(self.chosen_answer, "")
        
        if correct_skeleton not in chosen_name:
            self.errors.append(
                f"The chosen answer '{self.chosen_answer}' has an incorrect molecular skeleton. "
                f"The Diels-Alder product should be a '{correct_skeleton}', but the answer is a '{chosen_name}'."
            )

    def check_diels_alder_stereospecificity(self):
        """
        Checks if the chosen answer's stereochemistry is consistent with the
        stereochemistry of the starting materials.
        """
        # Constraint 1: The dienophile is 'cis-but-2-ene'.
        # This means the two methyl groups it provides (at C5 and C6) must be 'cis' to each other.
        
        # Extract R/S configurations for C5 and C6 from the chosen answer's name.
        chosen_name = self.options.get(self.chosen_answer)
        if not chosen_name or "4,5,6" not in chosen_name:
            return # Skeleton is wrong, no need to check stereochem.

        match = re.search(r'\(.*5(S|R),.*6(S|R)\)', chosen_name)
        if not match:
            self.errors.append(f"Could not parse C5/C6 stereochemistry from the name: {chosen_name}")
            return
            
        c5_config, c6_config = match.groups()
        
        # Determine the relative stereochemistry from the R/S labels.
        relative_stereochem = self.get_relative_stereochem_from_rs(c5_config, c6_config)
        
        if relative_stereochem != "cis":
            self.errors.append(
                f"Constraint violation: The starting dienophile is 'cis-but-2-ene', which requires the methyl groups at C5 and C6 of the product to be 'cis'."
            )
            self.errors.append(
                f"The chosen answer '{self.chosen_answer}' has a configuration of (5{c5_config}, 6{c6_config}), which corresponds to a 'trans' relationship. Therefore, the answer is incorrect."
            )

    def check_reasoning_consistency(self):
        """
        Checks if the provided reasoning is internally consistent and correct.
        """
        # The reasoning correctly states the rules for the endo product.
        # It then claims the absolute configuration is (1S, 4R, 5S, 6S).
        # Let's check if the reasoning's final R/S calculation is correct.
        if "The calculated absolute configuration for the endo product is **(1S, 4R, 5S, 6S)**" in self.reasoning:
            # The correct R/S assignment for the endo product is (1S, 4R, 5S, 6R).
            self.errors.append(
                "The reasoning contains a calculation error. The correct absolute configuration for the endo product described is (1S, 4R, 5S, 6R), which corresponds to option B, not C."
            )


    def run_checks(self):
        """Runs all checks and returns the final verdict."""
        self.check_product_skeleton()
        self.check_diels_alder_stereospecificity()
        self.check_reasoning_consistency()

        if not self.errors:
            return "Correct"
        else:
            # Use a set to get unique error messages and then join them.
            unique_errors = sorted(list(set(self.errors)))
            return "Incorrect. " + " ".join(unique_errors)

# --- Provided Data ---
llm_final_answer = "C"
llm_reasoning_text = """
Here is a step-by-step analysis of the chemical problem to determine the correct answer.
### **Step 1: Identification of Compound A**
...
### **Step 2: Identification of Compound B**
...
### **Step 3: Identification of Compound C**
...
### **Step 4: The Diels-Alder Reaction and Identification of Compound D**
...
*   **Product Skeleton**: The reaction forms a six-membered ring. The substituents from the diene (-OH and -CH₃) and the dienophile (two -CH₃ groups) result in a **4,5,6-trimethylcyclohex-2-enol** structure. This immediately eliminates options A and D, which have an incorrect 4,6,6-trimethyl substitution pattern. The correct answer must be either B or C.
### **Step 5: Stereochemical Analysis**
...
**Combining the Rules for the Endo Product:**
*   The {C1-OH, C4-CH₃} group is **cis**.
*   The {C5-CH₃, C6-CH₃} group is **cis**.
*   The {C1, C4} group is **trans** to the {C5, C6} group.
**Assigning Absolute Configuration (R/S) to the Endo Product:**
...
The calculated absolute configuration for the endo product is **(1S, 4R, 5S, 6S)**. This corresponds exactly to option C.
"""

# --- Execution ---
checker = ChemistryChecker(llm_final_answer, llm_reasoning_text)
result = checker.run_checks()
print(result)