import sys
from io import StringIO

class ChemistrySymmetryChecker:
    """
    This class checks the correctness of the answer to a multi-step synthesis and symmetry problem.
    """

    def __init__(self, question_details, provided_answer):
        self.question = question_details['question']
        self.options = question_details['options']
        self.final_answer_letter = provided_answer.strip('<>')
        self.final_answer_group = self.options.get(self.final_answer_letter)

    def get_product_info(self):
        """
        Analyzes the reaction and returns information about plausible products and their symmetries.
        """
        # Step 1: Toluene + HNO3/H2SO4 -> p-nitrotoluene (Product 1)
        # Step 2: p-nitrotoluene + MnO2/H2SO4 -> p-nitrobenzaldehyde (Product 2)
        # This pathway is chemically sound and necessary for the next step.

        # Step 3: p-nitrobenzaldehyde + acetone + NaOH -> Product 3
        # This is a Claisen-Schmidt condensation with two plausible outcomes.

        products = {
            "single_condensation": {
                "name": "(E)-4-(4-nitrophenyl)but-3-en-2-one",
                "point_group": "Cs",
                "comment": "This planar molecule is asymmetric end-to-end. Its only symmetry element besides identity is the molecular plane itself."
            },
            "double_condensation": {
                "name": "(1E,4E)-1,5-bis(4-nitrophenyl)penta-1,4-dien-3-one",
                "point_group": "C2v",
                "comment": "This molecule has a C2 axis passing through the C=O bond and two mirror planes containing this axis (one being the molecular plane). It does NOT have a center of inversion or perpendicular C2 axes due to the asymmetry of the C=O group."
            }
        }
        return products

    def check_correctness(self):
        """
        Checks if the provided answer is correct based on chemical principles.
        """
        products_info = self.get_product_info()
        double_condensation_product = products_info["double_condensation"]
        
        # The provided answer is 'C', which corresponds to 'd2h'.
        # The reasoning for this answer is based on the formation of the double condensation product.
        
        if self.final_answer_group is None:
            return f"Invalid answer format. The letter '{self.final_answer_letter}' does not correspond to any option."

        # Check if the provided answer's symmetry group matches the actual symmetry of the reasoned product.
        if self.final_answer_group == double_condensation_product['point_group']:
            # This case would be if the answer was C2v, which it is not.
            return "Correct"
        
        # The provided answer is D2h, but the actual point group of the molecule is C2v.
        if self.final_answer_group == "d2h":
            reason = (
                f"The provided answer 'd2h' is incorrect.\n"
                f"The reasoning correctly identifies the most plausible product as the double condensation product, "
                f"'(1E,4E)-1,5-bis(4-nitrophenyl)penta-1,4-dien-3-one'.\n"
                f"However, the symmetry analysis is flawed. This molecule's correct point group is C2v, not D2h.\n"
                f"A D2h point group requires, among other things, a center of inversion and three mutually perpendicular C2 axes. "
                f"The central carbonyl (C=O) group breaks this high symmetry; for example, there is no center of inversion at the central carbon atom because the oxygen atom on one side has no counterpart on the other.\n"
                f"Since the actual point group of the most likely product (C2v) is not an option, and the assigned point group (D2h) is incorrect, the answer is wrong."
            )
            return reason

        # Check against other possibilities
        if self.final_answer_group == "c2h":
             reason = (
                f"The provided answer 'c2h' is incorrect.\n"
                f"This answer is likely based on the double condensation product but with a flawed symmetry analysis. "
                f"The molecule '(1E,4E)-1,5-bis(4-nitrophenyl)penta-1,4-dien-3-one' does not have a center of inversion, which is a required element for the C2h point group. Its correct point group is C2v."
            )
             return reason

        if self.final_answer_group == "cs":
            reason = (
                f"The provided answer 'cs' is based on the single condensation product. While this product does have Cs symmetry, the provided final answer from the LLM is 'd2h', not 'cs'. Therefore, this check is for the provided answer 'd2h'."
            )
            # This is not the case for the final answer, but we include the logic.
            return "This check is not for the 'cs' answer."

        return f"The provided answer '{self.final_answer_group}' is incorrect for reasons related to the analysis above."


def run_check():
    """
    Runs the check for the given problem and answer.
    """
    # Based on the prompt, the options are:
    # A) c3, B) cs, C) d2h, D) c2h
    question_details = {
        "question": "Synthesis of Product 3 and its molecular symmetry group.",
        "options": {'A': 'c3', 'B': 'cs', 'C': 'd2h', 'D': 'c2h'}
    }

    # The final answer from the LLM to be checked
    provided_answer = "<<<C>>>"

    checker = ChemistrySymmetryChecker(question_details, provided_answer)
    result = checker.check_correctness()
    
    # Capture the output to display it
    old_stdout = sys.stdout
    sys.stdout = captured_output = StringIO()
    
    print(result)
    
    sys.stdout = old_stdout
    return captured_output.getvalue()

# Execute the check and print the result
output = run_check()
print(output)
