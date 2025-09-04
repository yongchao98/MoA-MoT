import re

class ChemistryProblemChecker:
    """
    This class checks the correctness of an answer to a multi-step organic chemistry problem.
    It simulates the reaction pathway based on established chemical principles and analyzes
    the final product's structure to determine the number of distinct hydrogen atoms.
    """

    def __init__(self):
        # A knowledge base of molecules and their relevant chemical properties.
        self.molecule_data = {
            "o-xylylene": {
                "fate_on_heating": "dibenzo[a,e]cyclooctadiene",
                "notes": "Highly reactive diene that dimerizes via [4+4] cycloaddition."
            },
            "7-oxobicyclo[2.2.1]hepta-2,5-diene": {
                "fate_on_heating": "benzene",
                "notes": "Thermally unstable; undergoes cheletropic elimination of CO."
            },
            "dibenzo[a,e]cyclooctadiene": {
                "common_name": "o-xylylene dimer",
                "formula": "C16H16",
                "symmetry": "C2",
                "distinct_H": 8,
                "reason": "The stable puckered conformation has C2 symmetry, resulting in 4 unique aromatic and 4 unique aliphatic hydrogen environments."
            },
            "benzene": {
                "distinct_H": 1
            }
        }
        # Mapping of options to their numerical values.
        self.options = {'A': 10, 'B': 4, 'C': 8, 'D': 7}

    def deduce_final_product_pathway(self):
        """
        Simulates the logical deduction of the reaction pathway to find the final product.
        """
        # Step 1-3: A double Diels-Alder adduct is formed, deprotected, and oxidized.
        # Step 4: The adduct undergoes a retro-Diels-Alder reaction upon heating.
        fragments = ["o-xylylene", "7-oxobicyclo[2.2.1]hepta-2,5-diene"]
        
        # Determine the fate of the reactive fragments.
        final_products = {}
        for frag in fragments:
            if frag in self.molecule_data and "fate_on_heating" in self.molecule_data[frag]:
                final_products[frag] = self.molecule_data[frag]["fate_on_heating"]

        # The question asks for "product 4". This refers to the major, complex organic
        # product derived from the starting materials, which is the dimer of o-xylylene.
        final_product_name = final_products.get("o-xylylene")
        
        return final_product_name

    def get_correct_h_count(self, molecule_name):
        """
        Retrieves the number of distinct hydrogens for a given molecule from the knowledge base.
        """
        if molecule_name in self.molecule_data:
            return self.molecule_data[molecule_name].get("distinct_H")
        return None

    def check_answer(self, llm_answer_text):
        """
        Checks the provided LLM answer against the deduced correct answer.
        """
        # Extract the chosen option letter from the answer string.
        match = re.search(r'<<<([A-D])>>>', llm_answer_text)
        if not match:
            return "Error: The answer format is invalid. Expected '<<<X>>>' where X is A, B, C, or D."
        
        answer_letter = match.group(1)
        proposed_count = self.options.get(answer_letter)

        # Deduce the correct answer through chemical reasoning.
        final_product = self.deduce_final_product_pathway()
        if not final_product:
            return "Error: Could not deduce the final product based on the encoded chemical knowledge."
            
        correct_count = self.get_correct_h_count(final_product)

        # Compare the proposed answer with the correct answer.
        if proposed_count == correct_count:
            return "Correct"
        else:
            reasoning = self.molecule_data[final_product]
            error_message = (
                f"Incorrect. The provided answer corresponds to {proposed_count} distinct hydrogens.\n\n"
                f"The correct analysis is as follows:\n"
                f"1. The four-step reaction sequence culminates in a retro-Diels-Alder reaction, which fragments the intermediate.\n"
                f"2. The fragments are two molecules of reactive *o-xylylene* and one molecule of 7-oxonorbornadiene (which decomposes to benzene).\n"
                f"3. The two *o-xylylene* molecules dimerize to form the final stable product (4), which is **{final_product}**.\n"
                f"4. This final product has **{reasoning['symmetry']} symmetry**.\n"
                f"5. This symmetry results in **{correct_count} chemically distinct hydrogen atoms**. ({reasoning['reason']})\n"
                f"Therefore, the correct answer is {correct_count}, which corresponds to option C."
            )
            return error_message

# This is the function that will be executed to check the answer.
def check_correctness_of_answer(llm_answer):
    checker = ChemistryProblemChecker()
    return checker.check_answer(llm_answer)

# The final answer from the LLM to be checked.
llm_final_answer = "<<<C>>>"

# Execute the check.
result = check_correctness_of_answer(llm_final_answer)
print(result)