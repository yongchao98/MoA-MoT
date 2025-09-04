import re

class ChemistryPuzzleChecker:
    """
    This class verifies the solution to the chemical puzzle by checking each logical step
    and calculation presented in the answer.
    """

    def __init__(self, llm_answer_text):
        """
        Initializes the checker with the necessary data and the LLM's answer.
        """
        self.llm_answer = llm_answer_text
        self.errors = []
        
        # Data from the problem statement and general chemistry knowledge
        self.options = {'A': 25, 'B': 31, 'C': 29, 'D': 35}
        self.atomic_data = {
            'H': {'Z': 1, 'mass_num': 1},
            'D': {'Z': 1, 'mass_num': 2},  # Deuterium, isotope of H
            'Li': {'Z': 3, 'mass_num': 7},
            'Al': {'Z': 13, 'mass_num': 27}
        }
        self.substance_properties = {
            'H2O': {'melting_point_K': 273.15},
            'D2O': {'melting_point_K': 276.97}
        }
        self.gas_properties = {
            'D2': {'protons': 2, 'neutrons': 2}
        }

    def check_reasoning(self):
        """
        Verifies the logical deductions made to identify the substances.
        """
        # 1. Check identification of Substance B (D2O)
        mp_b_approx = 277
        mp_d2o = self.substance_properties['D2O']['melting_point_K']
        if not abs(mp_b_approx - mp_d2o) < 1:
            self.errors.append(f"Reasoning Error: The melting point of B (~277 K) is not a close match for D2O ({mp_d2o} K).")

        # 2. Check identification of Gas W (D2)
        gas_w_protons = self.gas_properties['D2']['protons']
        gas_w_neutrons = self.gas_properties['D2']['neutrons']
        if gas_w_protons != gas_w_neutrons:
            self.errors.append("Reasoning Error: The identified gas D2 does not have an equal number of protons and neutrons, violating a key constraint.")

        # 3. The identification of Substance X as LiAlD4 is consistent with all clues.
        # This complex reasoning is accepted as correct based on the provided analysis.
        
    def perform_calculation(self):
        """
        Performs the calculation based on the identified Substance X (LiAlD4).
        """
        formula = "LiAlD4"
        elements_in_formula = ['Li', 'Al', 'D']

        # Identify lightest and heaviest elements by atomic number (Z)
        lightest_element_symbol = min(elements_in_formula, key=lambda el: self.atomic_data[el]['Z'])
        heaviest_element_symbol = max(elements_in_formula, key=lambda el: self.atomic_data[el]['Z'])

        if lightest_element_symbol != 'D' or heaviest_element_symbol != 'Al':
            self.errors.append(f"Calculation Error: Incorrectly identified lightest/heaviest elements. Expected D and Al, but got {lightest_element_symbol} and {heaviest_element_symbol}.")
            return None

        # Count atoms and sum masses based on the formula Li(1)Al(1)D(4)
        count_lightest = 4  # 4 atoms of D
        count_heaviest = 1  # 1 atom of Al

        mass_from_lightest = count_lightest * self.atomic_data[lightest_element_symbol]['mass_num']
        mass_from_heaviest = count_heaviest * self.atomic_data[heaviest_element_symbol]['mass_num']
        
        total_mass = mass_from_lightest + mass_from_heaviest
        
        # Verify the calculation steps mentioned in the LLM's answer
        if mass_from_lightest != 8 or mass_from_heaviest != 27 or total_mass != 35:
            self.errors.append(f"Calculation Error: The intermediate or final calculation is wrong. Expected 4*2 + 1*27 = 35.")
            return None
            
        return total_mass

    def check_final_answer(self):
        """
        Checks if the final selected option matches the calculated result.
        """
        # Extract the final answer choice from the LLM's response
        match = re.search(r'<<<([A-D])>>>', self.llm_answer)
        if not match:
            self.errors.append("Format Error: Could not find the final answer in the required '<<<X>>>' format.")
            return

        llm_choice_letter = match.group(1)
        
        # Get the correct calculated value
        calculated_value = self.perform_calculation()
        if calculated_value is None:
            return # Errors already logged

        # Find which option letter corresponds to the correct value
        correct_choice_letter = None
        for option, value in self.options.items():
            if value == calculated_value:
                correct_choice_letter = option
                break
        
        if correct_choice_letter is None:
            self.errors.append(f"Options Error: The calculated value {calculated_value} does not correspond to any of the given options.")
            return

        # Compare the LLM's choice with the correct choice
        if llm_choice_letter != correct_choice_letter:
            self.errors.append(f"Incorrect Answer: The final answer should be {calculated_value}, which is option {correct_choice_letter}, but the provided answer was {llm_choice_letter}.")

    def run_checks(self):
        """
        Runs all verification steps and returns the final result.
        """
        self.check_reasoning()
        self.check_final_answer()

        if self.errors:
            return "\n".join(self.errors)
        else:
            return "Correct"

# The provided answer from the LLM to be checked
llm_answer_to_check = """
Here is a step-by-step analysis to determine the correct answer based on the provided candidate responses.

**1. Analysis of the Chemical Puzzle**

A strong consensus exists among the detailed candidate answers (1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 15) regarding the identity of the substances involved. The reasoning is consistent and sound:
*   **Substance B:** Its melting point of ~277 K is a near-perfect match for heavy water (D₂O), which melts at 276.97 K.
*   **Gas W:** A molecule with an equal number of protons and neutrons, in the context of a reaction involving deuterium, is almost certainly deuterium gas (D₂), which has 2 protons and 2 neutrons.
*   **Substance X:** The clues point overwhelmingly to **Lithium Aluminum Deuteride (LiAlD₄)**.
    *   It contains a heavier isotope (Deuterium, D).
    *   Its non-deuterated analog, LiAlH₄, is a very common and powerful reducing agent in organic chemistry.
    *   It reacts violently with water (or in this case, heavy water, D₂O, which is liquid Y).
    *   The reaction LiAlD₄ + 4D₂O → LiOD + Al(OD)₃ + 4D₂ produces the correct gas (W=D₂) and a precipitate (G=Al(OD)₃).
    *   Heating the precipitate G (Al(OD)₃) releases substance B (D₂O).
    *   LiAlD₄ reduces a keto acid to a diol, a product with two oxygen atoms, matching the final clue.

**2. Analysis of the Calculation**

All candidate answers that perform the calculation follow the same correct logic:
*   **Substance X:** LiAlD₄
*   **Lightest Element:** Hydrogen (present as its isotope, Deuterium, D).
*   **Heaviest Element:** Aluminum (Al).
*   The question asks for the cumulative mass of *all* atoms of the lightest and heaviest elements.
*   **Mass from Lightest Element (Deuterium):** There are 4 atoms of D. Using the integer mass number, the total mass is 4 × 2 = 8.
*   **Mass from Heaviest Element (Aluminum):** There is 1 atom of Al. Using the integer mass number, the total mass is 1 × 27 = 27.
*   **Cumulative Sum:** 8 + 27 = **35**.

**3. Evaluation of Candidate Answers and Final Decision**

The calculated result is unequivocally 35. The provided options are:
A) 25
B) 31
C) 29
D) 35

*   A significant number of candidate answers (1, 2, 4, 9, 11) correctly calculate 35 but incorrectly select option C. Others (3, 6, 12, 15) calculate 35 but select option A. This indicates a simple mistake in matching the numerical result to the corresponding letter option.
*   Candidate answers 5 and 10 correctly perform the deduction, arrive at the value of 35, and correctly select option D.
*   Candidate answers 7, 8, 13, and 14 provide no reasoning, but 13 correctly selects D.

Based on the overwhelming consensus on the chemical identification and the mathematical calculation, the correct numerical answer is 35. This corresponds to option D. The errors in the other candidate answers are clerical, not logical.

<<<D>>>
"""

# Instantiate the checker and run the verification
checker = ChemistryPuzzleChecker(llm_answer_to_check)
result = checker.run_checks()
print(result)