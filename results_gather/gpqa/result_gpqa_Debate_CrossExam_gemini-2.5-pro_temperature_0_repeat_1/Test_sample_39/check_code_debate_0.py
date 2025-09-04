import re

def check_correctness_of_chemistry_answer(llm_response: str) -> str:
    """
    Checks the correctness of the LLM's answer to the chemistry jargon question.

    This function uses a rule-based approach, encoding common knowledge from a
    synthetic organic chemistry lab to validate the answer. It specifically checks
    the interpretation of the phrase "on top of each other".
    """

    # 1. Define the key phrase from the question and the options.
    key_phrase = "on top of each other"
    options = {
        "A": "Similar optical rotations",
        "B": "Non-covalent/van der Waals interactions",
        "C": "Similar boiling points",
        "D": "Similar polarities"
    }

    # 2. Create a knowledge base mapping lab jargon to scientific principles.
    # The 'match_strength' indicates how well the phrase describes the problem.
    jargon_interpretation = {
        "A": {
            "principle": "Similar optical rotations",
            "common_description": "This relates to stereochemistry (enantiomers). The phrase 'on top of each other' is not used to describe a failed separation, but rather the concept of superimposition.",
            "match_strength": 0.1
        },
        "B": {
            "principle": "Non-covalent/van der Waals interactions",
            "common_description": "This can cause aggregation or stickiness, but 'on top of each other' is not the standard term.",
            "match_strength": 0.2
        },
        "C": {
            "principle": "Similar boiling points",
            "common_description": "This refers to a failed distillation, resulting in co-distillation. The phrase is not a good fit.",
            "match_strength": 0.4
        },
        "D": {
            "principle": "Similar polarities",
            "common_description": "This is the classic description for a failed chromatography (TLC/column) separation, where unresolved spots literally appear 'on top of each other' on a TLC plate. This is a perfect fit.",
            "match_strength": 1.0
        }
    }

    # 3. Determine the most likely correct answer based on the knowledge base.
    best_option = max(jargon_interpretation, key=lambda k: jargon_interpretation[k]['match_strength'])

    # 4. Extract the LLM's final answer from the provided text.
    match = re.search(r'<<<([A-D])>>>', llm_response)
    if not match:
        return "Incorrect: The answer format is invalid. It should end with '<<<X>>>' where X is one of the options A, B, C, or D."

    llm_answer_key = match.group(1)

    # 5. Compare the LLM's answer with the determined correct answer.
    if llm_answer_key == best_option:
        return "Correct"
    else:
        reason = (
            f"Incorrect. The provided answer was <<<{llm_answer_key}>>>, but the correct answer is <<<{best_option}>>>.\n"
            f"Reasoning: The phrase '{key_phrase}' is a common idiom in organic chemistry labs that visually describes the result of a failed separation using chromatography (e.g., Thin-Layer Chromatography).\n"
            f"In chromatography, compounds are separated based on their polarity. If they have very similar polarities, they do not separate well, and their corresponding spots on a TLC plate overlap.\n"
            f"The chosen answer '{options[llm_answer_key]}' is a less likely interpretation of this specific phrase."
        )
        return reason

# The full response from the other LLM to be checked
llm_full_response = """
Proposer: The chemist's statement, "my compounds are on top of each other," is a common expression in a synthetic organic chemistry lab. It almost certainly refers to a failed separation attempt using chromatography, a primary technique for purifying compounds. The most common form, thin-layer chromatography (TLC) or column chromatography, separates compounds based on their polarity. If compounds have very similar polarities, they travel at the same rate through the stationary phase, resulting in poor separation. On a TLC plate, this would be visualized as spots that are not resolved and appear "on top of each other". Therefore, the chemist is most likely referring to the compounds having similar polarities. This aligns with option D.

Skeptic: While chromatography is a very common technique, it's not the only one. The chemist could have been attempting a distillation, and the compounds have similar boiling points (Option C), causing them to co-distill. This would also be a failed separation. Also, the term "on top of each other" could be interpreted literally in the context of stereochemistry, where one might try to superimpose mirror images (Option A). Furthermore, strong intermolecular forces (Option B) could cause compounds to aggregate, although this is less likely to be described in this specific way. The key is the common vernacular in a lab. While similar boiling points are a valid problem, the phrase "on top of each other" is most visually and idiomatically descriptive of a TLC plate where the spots have not separated.

Decision: The argument for chromatography is stronger due to the specific phrasing. The visual result of a TLC plate with unseparated compounds is literally spots on top of each other. This is a very frequent challenge and a common way to describe it. The search results confirm that the phrase "one spot on the other" or similar is used to describe poor separation in TLC due to close Rf values, which is a function of polarity. Therefore, the most probable answer is D.

<<<D>>>
"""

# Run the check
result = check_correctness_of_chemistry_answer(llm_full_response)
print(result)