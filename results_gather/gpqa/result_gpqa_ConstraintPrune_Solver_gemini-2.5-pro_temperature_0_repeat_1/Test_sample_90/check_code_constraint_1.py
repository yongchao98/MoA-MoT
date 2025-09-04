import sys
import io

class OrganicChemistryChecker:
    """
    This class checks the logic of a multi-step organic chemistry synthesis problem.
    It does not perform quantum calculations but encodes established reaction principles
    to verify the reasoning chain.
    """

    def __init__(self, question, llm_answer, llm_reasoning):
        self.question = question
        self.llm_answer = llm_answer
        self.llm_reasoning = llm_reasoning
        self.errors = []

    def check_rearrangement(self):
        """
        Checks if the reasoning correctly identifies that a rearrangement is necessary.
        """
        # Product 1 is 2-(hydroxy(phenyl)methyl)cyclohexan-1-one.
        # Direct fluorination of both functional groups would yield:
        # 1,1-difluoro-2-(fluoro(phenyl)methyl)cyclohexane.
        # The options A and C have the structure:
        # 2,2-difluoro-1-(fluoro(phenyl)methyl)cyclohexane.
        # These are constitutional isomers. A rearrangement must occur.
        if "rearrangement" not in self.llm_reasoning.lower():
            self.errors.append(
                "Constraint Violated: The reasoning fails to identify that a rearrangement is necessary. "
                "Direct fluorination of the aldol product leads to a different constitutional isomer "
                "than the one presented in options A and C."
            )
            return False
        return True

    def check_stereochemistry(self):
        """
        Checks the stereochemical argument from start to finish.
        This is the most critical part of the validation.
        """
        # Principle 1: Aldol Stereochemistry
        # Cyclohexanone + LDA (kinetic control) + Benzaldehyde -> anti-diastereomer is major.
        # Let's follow one enantiomer of the anti-product: (2R)-2-((R)-hydroxy(phenyl)methyl)cyclohexan-1-one.
        # So, the starting stereocenters are (Ring alpha-carbon: R, Benzylic alcohol: R).
        start_config = {'ring': 'R', 'benzylic': 'R'}

        # The LLM reasoning correctly assumes the anti-diastereomer.
        if "anti diastereomer" not in self.llm_reasoning:
            self.errors.append(
                "Reasoning Incomplete: The reasoning does not explicitly state the key assumption "
                "of forming the 'anti' diastereomer in the initial aldol reaction."
            )
            # We proceed assuming this is implied.

        # Principle 2: Fluorination Stereochemistry with DAST
        # DAST fluorination of a secondary alcohol typically proceeds with inversion of configuration (SN2).
        # The LLM reasoning correctly states this principle.
        if "inversion of configuration" not in self.llm_reasoning:
            self.errors.append(
                "Reasoning Flaw: The reasoning fails to mention the stereochemical outcome of the "
                "fluorination step (inversion), which is crucial for determining the final product."
            )
            return False

        # Apply the principles to predict the outcome:
        # Starting benzylic alcohol is 'R'.
        # After inversion, the benzylic fluoride should be 'S'.
        predicted_benzylic_config = 'S'

        # The stereocenter on the ring (migration origin) is assumed to retain its configuration.
        predicted_ring_config = start_config['ring'] # 'R'

        # Predicted final stereochemistry: (Ring: R, Benzylic: S)

        # Compare prediction with the options:
        # Option A: ((R)-((R)-...)-...) -> (Ring: R, Benzylic: R)
        # Option C: ((S)-((R)-...)-...) -> (Ring: R, Benzylic: S)
        
        # Our derived product (Ring: R, Benzylic: S) matches Option C.

        # Now, check the LLM's conclusion.
        if self.llm_answer == 'A':
            error_message = (
                "Incorrect: The model's reasoning is internally inconsistent. "
                "It correctly states two key principles: "
                "1. The reaction starts with the 'anti' aldol product (implying, for example, an R configuration at the benzylic alcohol center). "
                "2. The DAST fluorination of the alcohol proceeds with 'inversion of configuration'. "
                "Applying these principles correctly means the R alcohol becomes an S fluoride. "
                "This logical chain leads to option C, which has an S benzylic center. "
                "However, the model incorrectly concludes that this logic leads to option A, which has an R benzylic center. "
                "The conclusion contradicts the premises."
            )
            self.errors.append(error_message)
            return False
        
        if self.llm_answer == 'C' and not self.errors:
            # If the answer was C, the logic would be sound.
            return True
            
        return True


    def run_checks(self):
        """
        Runs all checks and returns the final result.
        """
        self.check_rearrangement()
        self.check_stereochemistry()

        if not self.errors:
            return "Correct"
        else:
            return "\n".join(self.errors)

def check_answer():
    """
    This function is the entry point for the checking process.
    """
    question = "cyclohexanone is treated with LDA at low temperature, followed by benzaldehyde. the reaction mixture is then acidified, forming product 1. 1 is then treated with an excess of diethylaminosulfur trifluoride, forming product 2. what is the structure of prodcut 2?"
    
    # This is the answer and reasoning provided by the other LLM.
    llm_answer = "A"
    llm_reasoning = """The model's reasoning is sound and correctly navigates a complex multi-step synthesis involving a non-obvious rearrangement.

*   **Constraint 1 & 2 (Initial Reactions):** The identification of the initial aldol product and the expected reactivity of excess DAST (fluorinating both the alcohol and the ketone) is correct.
*   **Constraint 3 (Rearrangement):** This is the most critical and difficult step. The model correctly deduced that a direct substitution/fluorination would lead to a different constitutional isomer than the ones presented in the options. The inference of a semi-pinacol type rearrangement is an excellent insight and is the key to solving the problem. This correctly prunes the candidate pool to A and C.
*   **Constraint 4 (Stereochemistry):** The stereochemical analysis is also correct. It correctly assumes the major *anti* diastereomer from the aldol reaction and, crucially, the standard **inversion of configuration** (S_N2 mechanism) for the fluorination of the secondary alcohol with DAST. This logic correctly singles out option A as the only one with the resulting stereochemistry.

The analysis is thorough, logical, and arrives at the correct conclusion based on established principles of organic chemistry. There are no errors in the reasoning."""

    checker = OrganicChemistryChecker(question, llm_answer, llm_reasoning)
    result = checker.run_checks()
    
    # The final output should be either "Correct" or a reason for being incorrect.
    print(result)

check_answer()