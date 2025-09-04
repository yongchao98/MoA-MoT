import re

def check_answer(question_data, final_answer):
    """
    Checks the correctness of the final answer by systematically verifying each step of the chemical reasoning.
    """
    
    # --- Step 1: Identify Compound A ---
    # The NMR data (6H triplet, 4H quartet) is a classic, albeit simplified, description for n-butane.
    # 6H (2xCH3) next to CH2 -> triplet. 4H (2xCH2) next to CH3 -> quartet (simplified).
    # This is the most plausible starting alkane for the subsequent reactions.
    compound_A = "n-butane"

    # --- Step 2: Identify Compound B ---
    # Monobromination of n-butane favors the more stable secondary radical, yielding 2-bromobutane.
    compound_B = "2-bromobutane"

    # --- Step 3: Identify Compound C ---
    # Elimination of 2-bromobutane with alc. KOH (Zaitsev's rule) yields but-2-ene.
    # The question states C has two geometrical isomers, which is true for but-2-ene (cis and trans).
    # If B were 1-bromobutane, the product would be but-1-ene, which has no geometric isomers.
    # This confirms the reaction pathway is correct.
    compound_C = "but-2-ene"
    dienophile = "cis-but-2-ene"

    # --- Step 4: Analyze the Diels-Alder Reaction ---
    
    # Define the expected product skeleton from the cycloaddition
    expected_skeleton = "4,5,6-trimethylcyclohex-2-enol"
    
    # Define the options from the question
    options = {
        "A": "(1S,4R,5S,6S)-4,5,6-trimethylcyclohex-2-enol",
        "B": "(1S,4R)-4,6,6-trimethylcyclohex-2-enol",
        "C": "(1S,4S)-4,6,6-trimethylcyclohex-2-enol",
        "D": "(1S,4R,5S,6R)-4,5,6-trimethylcyclohex-2-enol"
    }

    # Helper function to parse IUPAC names
    def parse_name(name):
        try:
            config_match = re.search(r'\((.*?)\)', name)
            skeleton_match = re.search(r'\)-(.*)', name)
            if not config_match or not skeleton_match:
                return None, None
            
            skeleton = skeleton_match.group(1)
            configs = {}
            for part in config_match.group(1).split(','):
                num = int(re.search(r'\d+', part).group())
                rs = re.search(r'[RS]', part).group()
                configs[num] = rs
            return skeleton, configs
        except (AttributeError, ValueError):
            return None, None

    # Helper function to check for cis relationship between adjacent stereocenters
    def are_adjacent_substituents_cis(config1, config2):
        # For adjacent stereocenters, (R,S) or (S,R) is cis. (R,R) or (S,S) is trans.
        return config1 != config2

    # --- Verification Logic ---
    
    # Check the provided final answer
    if final_answer not in options:
        return f"Incorrect. The final answer '{final_answer}' is not one of the valid options A, B, C, or D."

    answer_name = options[final_answer]
    answer_skeleton, answer_configs = parse_name(answer_name)

    # Constraint 1: Check the product skeleton
    if answer_skeleton != expected_skeleton:
        return (f"Incorrect. The answer '{final_answer}' has the skeleton '{answer_skeleton}', "
                f"but the Diels-Alder reaction between {dienophile} and (1E,3E)-penta-1,3-dien-1-ol "
                f"produces a '{expected_skeleton}' skeleton.")

    # Constraint 2: Check the stereochemistry from the dienophile
    # The dienophile is cis-but-2-ene, so the methyl groups it provides (at C5 and C6) must be cis.
    if 5 not in answer_configs or 6 not in answer_configs:
        return (f"Incorrect. The answer '{final_answer}' is missing stereochemical descriptors for C5 or C6, "
                f"which are essential for determining the product's stereochemistry.")

    c5_config = answer_configs[5]
    c6_config = answer_configs[6]

    if not are_adjacent_substituents_cis(c5_config, c6_config):
        return (f"Incorrect. The answer '{final_answer}' has a ({c5_config},{c6_config}) configuration at C5 and C6, "
                f"which corresponds to a TRANS relationship between the methyl groups. "
                f"The reaction uses a CIS-dienophile (cis-but-2-ene), which requires the substituents at C5 and C6 to be CIS.")

    # If all constraints are satisfied, the answer is correct.
    return "Correct"

# --- Run the check ---
# This data is extracted from the provided question and the LLM's final answer.
question_data = {
    "compound_A_nmr": "triplet at 0.9 ppm for 6H and a quartet at 1.3 ppm for 4H",
    "dienophile": "cis-isomer of C",
    "diene": "(1E,3E)-penta-1,3-dien-1-ol",
    "options": {
        "A": "(1S,4R,5S,6S)-4,5,6-trimethylcyclohex-2-enol",
        "B": "(1S,4R)-4,6,6-trimethylcyclohex-2-enol",
        "C": "(1S,4S)-4,6,6-trimethylcyclohex-2-enol",
        "D": "(1S,4R,5S,6R)-4,5,6-trimethylcyclohex-2-enol"
    }
}
# The final answer provided by the LLM to be checked.
llm_final_answer = "D"

# Get the result from the checking code.
result = check_answer(question_data, llm_final_answer)
print(result)