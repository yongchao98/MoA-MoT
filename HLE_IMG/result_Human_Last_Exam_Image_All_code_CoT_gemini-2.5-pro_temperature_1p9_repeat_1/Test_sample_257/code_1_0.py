import sys
import io

# A helper class to structure the reasoning process.
class NMR_Predictor:
    def __init__(self, start_material, reagents):
        self.start_material = start_material
        self.reagents = reagents
        self.product = "Compound 1"
        self.reasoning_steps = []

    def solve(self):
        """
        Executes the step-by-step analysis to determine the NMR properties.
        """
        self._analyze_reaction()
        self._analyze_nmr()
        self._print_results()

    def _analyze_reaction(self):
        """
        Step 1: Predict the structure of the product, Compound 1.
        """
        self.reasoning_steps.append("Step 1: Predicting the structure of Compound 1")
        self.reasoning_steps.append(
            f"The reaction of {self.start_material} with {self.reagents} is an electrophilic aromatic substitution."
        )
        self.reasoning_steps.append(
            "Concentrated H2SO4 for 16 hours causes sulfonation, adding -SO3H groups."
        )
        self.reasoning_steps.append(
            "The product is described as water-soluble, which is consistent with the presence of polar -SO3H groups."
        )
        self.reasoning_steps.append(
            "The most likely sites for substitution are the positions para to the activating N-propyl groups on the outer rings."
        )
        self.reasoning_steps.append(
            "This results in a symmetrical, disulfonated product."
        )

    def _analyze_nmr(self):
        """
        Step 2: Analyze the 1H NMR of Compound 1 to find the most deshielded proton.
        """
        self.reasoning_steps.append("\nStep 2: Identifying the most deshielded proton in Compound 1")
        self.reasoning_steps.append(
            "The most deshielded proton is typically in the most electron-poor region of a molecule."
        )
        self.reasoning_steps.append(
            "In this cationic, polycyclic aromatic system, the single proton on the central ring is in the most electron-deficient environment, similar to the H-9 proton of an acridinium salt."
        )
        self.reasoning_steps.append(
            "Therefore, this proton will appear at the highest chemical shift (most deshielded)."
        )

        """
        Step 3: Determine the splitting pattern and integration for this proton.
        """
        self.reasoning_steps.append("\nStep 3: Determining the splitting pattern and integration")
        # Integration
        self.integration_val = 1
        self.reasoning_steps.append(
            f"Integration: There is only one such proton in the structure. Thus, its signal integrates to {self.integration_val}H."
        )
        # Splitting
        self.adjacent_protons = 0
        self.reasoning_steps.append(
            f"Splitting: The splitting pattern is determined by the n+1 rule, where 'n' is the number of adjacent protons."
        )
        self.reasoning_steps.append(
            f"This proton has no protons on adjacent carbons, so n = {self.adjacent_protons}."
        )
        self.splitting_result = self.adjacent_protons + 1
        self.reasoning_steps.append(
            f"Calculation: n + 1 = {self.adjacent_protons} + 1 = {self.splitting_result}."
        )
        self.pattern = "singlet"
        self.reasoning_steps.append(
            f"A signal with {self.splitting_result} peak is called a {self.pattern}."
        )

    def _print_results(self):
        """
        Prints the logical steps and the final conclusion.
        """
        print("Here is the step-by-step reasoning to solve the problem:\n")
        for step in self.reasoning_steps:
            print(step)

        print("\n--------------------------------------------------------------")
        print("Final Answer:")
        print(
            "The highest deshielded proton peak in the 1H NMR of Compound 1 is a:"
        )
        print(f"\nSplitting Pattern: {self.pattern.capitalize()}")
        print(f"Integration: {self.integration_val}H")
        print("--------------------------------------------------------------")


# Main execution block
if __name__ == "__main__":
    problem_solver = NMR_Predictor(
        start_material="di-propyl diazaoxatriangulenium (Pr-DAOTA)",
        reagents="Conc. H2SO4, 16h"
    )
    problem_solver.solve()