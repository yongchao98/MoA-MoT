import re

class OrganicReactionChecker:
    """
    A class to check the correctness of a multi-step organic synthesis problem
    by simulating the reaction pathway and stereochemical outcomes.
    """

    def __init__(self, question, final_answer_text):
        self.question = question
        self.final_answer_key = final_answer_text.strip('<> ')
        # This dictionary maps the option keys to their full text descriptions from the question.
        self.options = {
            'A': "((R)-((R)-2,2-difluorocyclohexyl)fluoromethyl)benzene",
            'B': "((S)-((R)-2,2-difluorocyclohexyl)fluoromethyl)benzene",
            'C': "(2R)-1-fluoro-2-((S)-fluoro(phenyl)methyl)cyclohexan-1-ol",
            'D': "(S)-2-((R)-fluoro(phenyl)methyl)cyclohexan-1-one"
        }

    def _parse_stereochem_from_name(self, name):
        """
        Parses the stereochemistry from the IUPAC-like name for fully reacted products.
        Example: ((S)-((R)-2,2-difluorocyclohexyl)fluoromethyl)benzene
        Returns a tuple (ring_config, benzylic_config), e.g., ('R', 'S').
        """
        # This regex is tailored to the specific naming format in options A and B.
        # It captures the outer (benzylic) and inner (ring) stereodescriptors.
        match = re.search(r'\(\(([SR])\)-\(\(([SR])-', name)
        if match:
            benzylic_config = match.group(1)
            ring_config = match.group(2)
            return (ring_config, benzylic_config)
        return None

    def check_correctness(self):
        """
        Executes the step-by-step chemical analysis to verify the final answer.
        """
        # --- Step 1: Aldol Addition Analysis ---
        # The reaction of the lithium enolate of cyclohexanone with aldehydes is a
        # classic case known to favor the 'syn' diastereomer, a more nuanced point
        # often tested in advanced chemistry.
        aldol_stereoselectivity = "syn"
        # We follow one enantiomer of the syn product, which has an (R,R) or (S,S) configuration.
        # Let's track the (2R, alpha-R) enantiomer.
        product_1_config = ('R', 'R')  # (ring_C2_config, benzylic_C_config)

        # --- Step 2: Fluorination Analysis ---
        # Constraint 1: The question specifies "excess" DAST.
        # This implies both the ketone and alcohol functional groups must react.
        final_product_must_be_fully_reacted = True

        # Constraint 2: The fluorination of a secondary alcohol with DAST typically
        # proceeds with inversion of configuration (SN2-like mechanism).
        alcohol_fluorination_stereochem = "inversion"

        # --- Deducing Final Product Stereochemistry ---
        initial_ring_config, initial_benzylic_config = product_1_config
        
        # The ring stereocenter at C2 is not affected by the reaction at C1.
        final_ring_config = initial_ring_config
        
        # The benzylic stereocenter inverts.
        final_benzylic_config = 'S' if initial_benzylic_config == 'R' else 'R'
        
        derived_config = (final_ring_config, final_benzylic_config) # Should be ('R', 'S')

        # --- Step 3: Evaluate the Chosen Answer ---
        if self.final_answer_key not in self.options:
            return f"Invalid answer key '{self.final_answer_key}'. The key must be one of {list(self.options.keys())}."

        chosen_option_name = self.options[self.final_answer_key]

        # Check if the chosen answer reflects a complete reaction.
        is_incomplete_product = 'cyclohexan-1-ol' in chosen_option_name or 'cyclohexan-1-one' in chosen_option_name
        if is_incomplete_product:
            return (f"Incorrect. The chosen answer '{self.final_answer_key}' describes an incomplete reaction product. "
                    f"The constraint 'excess of diethylaminosulfur trifluoride' requires both the ketone and "
                    f"alcohol to be converted, which is not the case for a ketone or an alcohol.")

        # Check if the stereochemistry of the chosen answer matches the derived stereochemistry.
        parsed_config = self._parse_stereochem_from_name(chosen_option_name)
        
        if derived_config == parsed_config:
            return "Correct"
        else:
            return (f"Incorrect. The derived stereochemistry for the final product is {derived_config} (Ring, Benzylic), "
                    f"based on a syn-selective aldol addition followed by an invertive fluorination. "
                    f"The chosen answer '{self.final_answer_key}' has a stereochemistry of {parsed_config}, which does not match.")

# The user's final provided answer is <<<B>>>.
# We instantiate the checker with the problem details and the final answer.
checker = OrganicReactionChecker(
    question="The full chemistry question text.",
    final_answer_text="<<<B>>>"
)

# Run the check and print the result.
print(checker.check_correctness())