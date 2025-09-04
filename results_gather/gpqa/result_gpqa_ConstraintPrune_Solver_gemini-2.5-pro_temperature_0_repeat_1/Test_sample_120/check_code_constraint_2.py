import re

def check_epoxide_opening_answer():
    """
    Checks the correctness of the answer for the reaction of
    (1R,3R,4R,6S)-1,3,4-trimethyl-7-oxabicyclo[4.1.0]heptane with Me2CuLi.
    """
    # --- Problem Definition ---
    # Reactant: (1R,3R,4R,6S)-1,3,4-trimethyl-7-oxabicyclo[4.1.0]heptane
    # This structure has an epoxide between C1 and C6 of a cyclohexane ring.
    reactant = {
        'epoxide_carbons': ('C1', 'C6'),
        'substituents': {
            'C1': ['Me'],  # C1 is tertiary and has a methyl group
            'C6': ['H']   # C6 is secondary
        },
        'initial_stereochemistry': {
            'C1': 'R',
            'C3': 'R',
            'C4': 'R',
            'C6': 'S'
        }
    }
    # The provided answer from the other LLM
    llm_answer_choice = 'C'
    llm_answer_name = '(1R,2R,4R,5R)-1,2,4,5-tetramethylcyclohexan-1-ol'

    # --- Step 1: Verify Regioselectivity ---
    # The nucleophile (Me-) from the organocuprate attacks the less sterically hindered carbon.
    # C1 is a tertiary carbon with a methyl group. C6 is a secondary carbon.
    # Therefore, C6 is less hindered.
    attack_site = 'C6'
    
    # The LLM's reasoning also states attack at C6. This is consistent.

    # --- Step 2: Verify Stereoselectivity ---
    # The Sₙ2 attack causes inversion of configuration at the attack site.
    original_config = reactant['initial_stereochemistry'][attack_site]
    if original_config == 'S':
        inverted_config = 'R'
    elif original_config == 'R':
        inverted_config = 'S'
    else:
        return f"Invalid initial stereochemistry '{original_config}' for attack site {attack_site}."

    # The LLM's reasoning states the configuration at C6 inverts from S to R. This is consistent.
    # Expected config at the attacked carbon is 'R'.

    # --- Step 3: Determine Product Structure and IUPAC Nomenclature ---
    # The reaction opens the epoxide to form a cyclohexanol.
    # - The -OH group is on the original C1.
    # - A new methyl group is on the original C6.
    # - Original methyl groups are on C1, C3, and C4.
    # For IUPAC naming:
    # - The carbon with the -OH group (original C1) becomes the new C1.
    # - The ring is numbered towards the new methyl group (original C6) to give lowest locants.
    # This results in the following mapping from original carbons to new IUPAC numbers:
    # New C1 <- Original C1
    # New C2 <- Original C6
    # New C3 <- Original C5
    # New C4 <- Original C4
    # New C5 <- Original C3
    # New C6 <- Original C2
    
    # The four methyl groups are located at:
    # - Original C1 -> New C1
    # - Original C3 -> New C5
    # - Original C4 -> New C4
    # - The new methyl group at Original C6 -> New C2
    derived_locants = sorted([1, 2, 4, 5])

    # --- Step 4: Determine Final Stereochemistry based on the new numbering ---
    derived_stereochemistry = {}
    # New C1 (from Original C1): Configuration is retained.
    derived_stereochemistry['1'] = reactant['initial_stereochemistry']['C1']
    # New C2 (from Original C6): Configuration is inverted.
    derived_stereochemistry['2'] = inverted_config
    # New C4 (from Original C4): Configuration is retained.
    derived_stereochemistry['4'] = reactant['initial_stereochemistry']['C4']
    # New C5 (from Original C3): Configuration is retained.
    derived_stereochemistry['5'] = reactant['initial_stereochemistry']['C3']

    # --- Step 5: Parse and Compare with the LLM's Answer ---
    try:
        # Parse the locants from the name, e.g., "1,2,4,5"
        locants_match = re.search(r'\)-([\d,]+)-', llm_answer_name)
        llm_locants = sorted([int(l) for l in locants_match.group(1).split(',')])
        
        # Parse the stereodescriptors, e.g., "(1R,2R,4R,5R)"
        stereo_match = re.search(r'\((.*?)\)', llm_answer_name)
        llm_stereo_pairs = stereo_match.group(1).split(',')
        llm_stereochemistry = {pair[0]: pair[1] for pair in llm_stereo_pairs}
    except (AttributeError, IndexError, ValueError):
        return f"Failed to parse the provided answer name: '{llm_answer_name}'. The format is unexpected."

    # Check locants
    if derived_locants != llm_locants:
        return f"Constraint check failed: Incorrect substituent locants. The reaction should produce a 1,2,4,5-tetramethylcyclohexanol. Derived locants are {derived_locants}, but the answer implies {llm_locants}."

    # Check stereochemistry
    if derived_stereochemistry != llm_stereochemistry:
        return f"Constraint check failed: Incorrect stereochemistry. Derived stereochemistry is {derived_stereochemistry}, but the answer specifies {llm_stereochemistry}."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_epoxide_opening_answer()
print(result)