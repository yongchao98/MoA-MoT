import re

def check_answer_correctness():
    """
    Checks the correctness of the LLM's answer for the multi-step synthesis problem.
    It verifies two key stereochemical constraints for the major product:
    1. The two ester groups must be 'cis' due to the use of maleic anhydride.
    2. The two main bridges (ethano-dicarboxylate and methano) should be 'anti'
       to each other, as the second Diels-Alder attack occurs from the sterically
       less hindered face.
    """
    
    # Data from the question
    problem_data = {
        'A': {
            "name": "dimethyl (1R,4S,4aR,4bR,5R,8S,8aS,8bS,10S,11R)-1,4,4a,4b,5,8,8a,8b-octahydro-1,4-ethano-5,8-methanobiphenylene-10,11-dicarboxylate",
        },
        'B': {
            "name": "dimethyl (1S,4R,4aR,4bR,5S,8R,8aS,8bS,10S,11R)-1,4,4a,4b,5,8,8a,8b-octahydro-1,4-ethano-5,8-methanobiphenylene-10,11-dicarboxylate",
        },
        'C': {
            "name": "dimethyl (1R,4S,4aS,4bR,5R,8S,8aS,8bR,10R,11R)-1,4,4a,4b,5,8,8a,8b-octahydro-1,4-ethano-5,8-methanobiphenylene-10,11-dicarboxylate",
        },
        'D': {
            "name": "dimethyl (1R,4S,4aS,4bR,5S,8R,8aS,8bR,10S,11R)-1,4,4a,4b,5,8,8a,8b-octahydro-1,4-ethano-5,8-methanobiphenylene-10,11-dicarboxylate",
        }
    }
    
    # The final answer provided by the LLM to be checked
    llm_answer_key = 'A'
    
    selected_option = problem_data.get(llm_answer_key)
    if not selected_option:
        return f"Invalid answer key '{llm_answer_key}'. Valid options are A, B, C, D."

    name = selected_option["name"]

    # --- Constraint 1: Ester groups must be 'cis' ---
    ester_match = re.search(r'(10[RS]),(11[RS])', name)
    if not ester_match:
        return "Could not determine ester stereochemistry from the IUPAC name."
    
    ester_c10 = ester_match.group(1)[-1]
    ester_c11 = ester_match.group(2)[-1]
    
    if ester_c10 == ester_c11:
        ester_stereo = "trans"
    else:
        ester_stereo = "cis"

    if ester_stereo != "cis":
        return (f"Incorrect. The reaction starts with maleic anhydride, a cis-dienophile. "
                f"This requires the two ester groups (at C10 and C11) in the final product to be 'cis' (i.e., R,S or S,R). "
                f"The selected answer '{llm_answer_key}' has a '{ester_c10},{ester_c11}' configuration, which is '{ester_stereo}'.")

    # --- Constraint 2: Bridge orientation must be 'anti' for the major product ---
    # 'anti' adduct has trans-fusions at the central cyclobutane ring (4a,4b and 8a,8b).
    # 'syn' adduct has cis-fusions.
    # Cis fusion: R,R or S,S. Trans fusion: R,S or S,R.
    
    fusion_match = re.search(r'(4a[RS]),(4b[RS]),.*? (8a[RS]),(8b[RS])', name)
    if not fusion_match:
        return "Could not determine bridge orientation from the IUPAC name's fusion stereochemistry."

    fusion_4a = fusion_match.group(1)[-1]
    fusion_4b = fusion_match.group(2)[-1]
    fusion_8a = fusion_match.group(3)[-1]
    fusion_8b = fusion_match.group(4)[-1]

    side1_is_cis = (fusion_4a == fusion_4b)
    side2_is_cis = (fusion_8a == fusion_8b)

    if side1_is_cis and side2_is_cis:
        bridge_orientation = "syn"
    elif not side1_is_cis and not side2_is_cis:
        bridge_orientation = "anti"
    else:
        bridge_orientation = "mixed/unclear"

    if bridge_orientation != "anti":
        return (f"Incorrect. The major product is formed via the most sterically favorable pathway. "
                f"This involves the cyclopentadiene adding to the face *opposite* the existing ethano-dicarboxylate bridge, "
                f"resulting in an 'anti' configuration of the two bridges. "
                f"Analysis of the IUPAC name for answer '{llm_answer_key}' shows it has a '{bridge_orientation}' configuration "
                f"(based on its {fusion_match.group(1)}-{fusion_match.group(2)} and {fusion_match.group(3)}-{fusion_match.group(4)} fusions). "
                f"This contradicts the requirement for the major product.")

    return "Correct"

# Execute the check
result = check_answer_correctness()
print(result)