import re

class OrganicSynthesisChecker:
    """
    A class to simulate the described organic synthesis step-by-step
    and verify the final product against a given answer.
    """
    def __init__(self):
        # Initial state representing (R)-(+)-Limonene
        self.molecule = {
            'name': '(R)-(+)-Limonene',
            'stereocenters': {'C4': 'R'},
            'substituents': {'C1': 'methyl', 'C4': 'isopropenyl'},
            'double_bonds': ['endocyclic', 'exocyclic'],
            'errors': []
        }
        self.log = []

    def _log_step(self, step_name, reasoning, expected_stereochem):
        """Logs the outcome of a reaction step."""
        is_correct = self.molecule['stereocenters'] == expected_stereochem
        self.log.append({
            'step': step_name,
            'reasoning': reasoning,
            'outcome_stereochem': self.molecule['stereocenters'],
            'is_correct': is_correct
        })
        if not is_correct:
            self.molecule['errors'].append(
                f"Error in step '{step_name}': "
                f"Expected stereochemistry {expected_stereochem}, "
                f"but derived {self.molecule['stereocenters']}."
            )

    def step1_hydrogenation(self):
        """Simulates selective hydrogenation."""
        reasoning = "H2/PdC selectively reduces the less substituted exocyclic double bond."
        self.molecule['substituents']['C4'] = 'isopropyl'
        self.molecule['double_bonds'] = ['endocyclic']
        self.molecule['name'] = '(R)-4-isopropyl-1-methylcyclohex-1-ene'
        self._log_step("1. Hydrogenation", reasoning, {'C4': 'R'})

    def step2_epoxidation(self):
        """Simulates diastereoselective epoxidation."""
        reasoning = "m-CPBA attacks anti to the bulky C4(R) group, forming the (1S, 2R, 4R)-epoxide."
        self.molecule['double_bonds'] = []
        self.molecule['stereocenters'].update({'C1': 'S', 'C2': 'R'})
        self.molecule['name'] = '(1S, 2R, 4R)-epoxide'
        self._log_step("2. Epoxidation", reasoning, {'C1': 'S', 'C2': 'R', 'C4': 'R'})

    def step3_epoxide_opening(self):
        """Simulates nucleophilic epoxide opening."""
        reasoning = "Basic NaOMe attacks the less hindered C2 via SN2, causing inversion (R -> S)."
        if self.molecule['stereocenters'].get('C2') == 'R':
            self.molecule['stereocenters']['C2'] = 'S'
        self.molecule['name'] = '(1S, 2S, 4R)-alcohol'
        self._log_step("3. Epoxide Opening", reasoning, {'C1': 'S', 'C2': 'S', 'C4': 'R'})

    def step4_esterification(self):
        """Simulates Steglich esterification."""
        reasoning = "Steglich esterification occurs with retention of all stereocenters."
        self.molecule['name'] = '(1S, 2S, 4R)-4-isopropyl-2-methoxy-1-methylcyclohexyl propionate'
        self._log_step("4. Esterification", reasoning, {'C1': 'S', 'C2': 'S', 'C4': 'R'})

    def run_synthesis(self):
        """Executes the full reaction sequence."""
        self.step1_hydrogenation()
        self.step2_epoxidation()
        self.step3_epoxide_opening()
        self.step4_esterification()
        return self.molecule

def check_correctness_of_answer():
    """
    Checks if the provided LLM answer is correct by simulating the synthesis
    and comparing the result to the chosen option.
    """
    # The final answer from the LLM to be checked
    llm_answer_choice = "D"

    # The options as stated in the question
    options = {
        "A": "1-methoxy-2-((S)-4-methylcyclohex-3-en-1-yl)propan-2-yl propionate",
        "B": "(1S,2R,4R)-4-isopropyl-2-methoxy-1-methylcyclohexyl propionate",
        "C": "(1S,2S,5R)-5-isopropyl-2-methoxy-2-methylcyclohexyl propionate",
        "D": "(1S,2S,4R)-4-isopropyl-2-methoxy-1-methylcyclohexyl propionate"
    }

    # 1. Simulate the synthesis to derive the correct product
    simulator = OrganicSynthesisChecker()
    derived_product = simulator.run_synthesis()
    derived_product_name = derived_product['name']

    # Check for any errors during the simulation
    if derived_product['errors']:
        return f"Incorrect. The simulation failed. Error: {derived_product['errors'][0]}"

    # 2. Get the product name corresponding to the LLM's answer
    chosen_option_name = options.get(llm_answer_choice)
    if not chosen_option_name:
        return f"Incorrect. The answer choice '{llm_answer_choice}' is not a valid option."

    # 3. Compare the derived product with the chosen option
    # We will check both the full name and the key stereochemical descriptor
    derived_stereochem = re.search(r'\((.*?)\)', derived_product_name).group(1)
    chosen_stereochem_match = re.search(r'\((.*?)\)', chosen_option_name)
    
    if not chosen_stereochem_match:
        return (f"Incorrect. The chosen option '{llm_answer_choice}' is structurally wrong and lacks a "
                f"standard stereochemical descriptor to compare. The correct product is {derived_product_name}.")

    chosen_stereochem = chosen_stereochem_match.group(1)

    if derived_product_name == chosen_option_name:
        return "Correct"
    else:
        error_reasons = []
        if derived_stereochem != chosen_stereochem:
            error_reasons.append(f"the stereochemistry is incorrect (derived {derived_stereochem}, but option {llm_answer_choice} has {chosen_stereochem})")
        
        # Add other structural checks if needed
        if "4-isopropyl" not in chosen_option_name:
             error_reasons.append("the isopropyl group is at the wrong position (should be C4)")
        
        if not error_reasons:
             error_reasons.append("the overall structure or naming is incorrect")

        return (f"Incorrect. The chosen option '{llm_answer_choice}' does not match the derived product because "
                f"{' and '.join(error_reasons)}. The correct product is: {derived_product_name}.")

# Run the check and print the result
print(check_correctness_of_answer())