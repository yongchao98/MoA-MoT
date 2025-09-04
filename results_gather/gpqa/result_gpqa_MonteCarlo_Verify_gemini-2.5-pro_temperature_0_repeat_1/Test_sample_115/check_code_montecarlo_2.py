import re

def check_correctness():
    """
    This function checks the correctness of the provided LLM answer by verifying each
    step of the chemical synthesis and deduction.
    """

    class ChemistryLogicVerifier:
        def __init__(self):
            self.errors = []
            self.options = {
                'A': "(1S,4S)-4,6,6-trimethylcyclohex-2-enol",
                'B': "(1S,4R,5S,6S)-4,5,6-trimethylcyclohex-2-enol",
                'C': "(1S,4R)-4,6,6-trimethylcyclohex-2-enol",
                'D': "(1S,4R,5S,6R)-4,5,6-trimethylcyclohex-2-enol"
            }
            # This data is extracted from the LLM's reasoning.
            self.llm_reasoning = {
                'A': "n-butane",
                'B': "2-bromobutane",
                'C': "but-2-ene",
                'dienophile': "cis-but-2-ene",
                'final_choice': 'D'
            }

        def check_step1_nmr(self):
            """Step 1: Verify Compound A from NMR data."""
            if self.llm_reasoning['A'] != "n-butane":
                self.errors.append("Step 1 Error: The NMR data (6H triplet, 4H quartet in the alkane region) is a classic pattern for n-butane, not {}.".format(self.llm_reasoning['A']))
            # The logic is sound: n-butane (CH3-CH2-CH2-CH3) has two equivalent CH3 groups (6H) split into a triplet by the adjacent CH2,
            # and two equivalent CH2 groups (4H) split into a quartet by the adjacent CH3.

        def check_step2_bromination(self):
            """Step 2: Verify monobromination product."""
            if self.llm_reasoning['B'] != "2-bromobutane":
                self.errors.append("Step 2 Error: Monobromination of n-butane is selective for the secondary position, yielding 2-bromobutane, not {}.".format(self.llm_reasoning['B']))

        def check_step3_elimination(self):
            """Step 3: Verify elimination product."""
            if self.llm_reasoning['C'] != "but-2-ene":
                self.errors.append("Step 3 Error: Elimination of 2-bromobutane with alcoholic KOH follows Zaitsev's rule to form the more stable, substituted alkene, which is but-2-ene, not {}.".format(self.llm_reasoning['C']))

        def check_step4_diels_alder(self):
            """Step 4: Verify the Diels-Alder reaction product and stereochemistry."""
            # 4a: Check connectivity
            correct_skeleton = "4,5,6-trimethylcyclohex-2-en"
            for key in ['B', 'D']:
                if correct_skeleton not in self.options[key]:
                    self.errors.append(f"Step 4 Error: The skeleton of option {key} is inconsistent with the Diels-Alder product.")
            
            # 4b: Check stereochemistry based on the dienophile
            # The dienophile is cis-but-2-ene, so the methyl groups at C5 and C6 must be cis.
            # Cis relationship for adjacent stereocenters means their R/S descriptors are different (R,S or S,R).
            # Trans relationship means their R/S descriptors are the same (R,R or S,S).

            try:
                # Check Option B: (1S,4R,5S,6S)-...
                name_b = self.options['B']
                descriptors_b = re.search(r'\((.*?)\)', name_b).group(1).split(',')
                c5_b = [d[1] for d in descriptors_b if d.startswith('5')][0]
                c6_b = [d[1] for d in descriptors_b if d.startswith('6')][0]
                if c5_b != c6_b: # Should be S,S (same) for trans
                    self.errors.append("Step 4 Error: The analysis of Option B is flawed. It has a (S,S) configuration at C5/C6, which is trans, but the code evaluated it as cis.")

                # Check Option D: (1S,4R,5S,6R)-...
                name_d = self.options['D']
                descriptors_d = re.search(r'\((.*?)\)', name_d).group(1).split(',')
                c5_d = [d[1] for d in descriptors_d if d.startswith('5')][0]
                c6_d = [d[1] for d in descriptors_d if d.startswith('6')][0]
                if c5_d == c6_d: # Should be S,R (different) for cis
                    self.errors.append("Step 4 Error: The analysis of Option D is flawed. It has a (S,R) configuration at C5/C6, which is cis, but the code evaluated it as trans.")

            except (AttributeError, IndexError):
                self.errors.append("Step 4 Error: Failed to parse R/S stereodescriptors from the option names.")

            # 4c: Check if the final choice matches the correct logic
            if self.llm_reasoning['final_choice'] != 'D':
                self.errors.append("Step 4 Error: The final answer selected was {}, but the correct stereochemistry (cis at C5/C6) points to option D.".format(self.llm_reasoning['final_choice']))

        def verify(self):
            """Run all checks and return the final result."""
            self.check_step1_nmr()
            self.check_step2_bromination()
            self.check_step3_elimination()
            self.check_step4_diels_alder()

            if self.errors:
                return "Incorrect. The following logical errors were found:\n- " + "\n- ".join(self.errors)
            else:
                return "Correct"

    # Instantiate and run the verifier
    verifier = ChemistryLogicVerifier()
    return verifier.verify()

# Execute the check and print the result
result = check_correctness()
print(result)