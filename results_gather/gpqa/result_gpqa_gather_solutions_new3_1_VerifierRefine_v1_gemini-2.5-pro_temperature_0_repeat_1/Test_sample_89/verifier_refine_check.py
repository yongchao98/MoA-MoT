import re

def check_correctness_of_answer():
    """
    This function checks the correctness of the final answer by comparing its
    IUPAC name against the properties of the product derived from a step-by-step
    chemical analysis.
    """
    # The final answer provided in the prompt is D, which corresponds to the name:
    llm_provided_name = "3,4-dimethyl-5,6-dioxooctanoic acid"

    # --- Ground Truth Properties from Chemical Analysis ---
    # The final product is HOOC-CH2-CH(CH3)-CH(CH3)-C(=O)-C(=O)-CH2-CH3
    correct_properties = {
        "parent_stem": "oct",
        "parent_length": 8,
        "principal_group": "oic acid",
        "substituents": {
            "methyl": sorted([3, 4]),
            "oxo": sorted([5, 6])
        },
        "total_carbons": 10  # 8 from parent + 2 from methyls
    }

    # --- Parse the LLM's Provided Name ---
    try:
        # 1. Parse parent chain and principal group
        match = re.search(r'(\w+)an(oic acid|al|one|trione)$', llm_provided_name)
        if not match:
            return f"Could not parse the parent chain and principal group from the name: '{llm_provided_name}'"

        parent_map = {"meth": 1, "eth": 2, "prop": 3, "but": 4, "pent": 5, "hex": 6, "hept": 7, "oct": 8, "non": 9}
        parsed_stem = match.group(1)
        parsed_group = match.group(2)
        parsed_length = parent_map.get(parsed_stem)

        # 2. Check parent chain and principal group
        if parsed_stem != correct_properties["parent_stem"]:
            return f"Incorrect parent chain. Expected '{correct_properties['parent_stem']}' (length {correct_properties['parent_length']}), but found '{parsed_stem}' (length {parsed_length})."
        if parsed_group != correct_properties["principal_group"]:
            return f"Incorrect principal functional group. Expected '{correct_properties['principal_group']}', but found '{parsed_group}'."

        # 3. Parse substituents
        substituent_part = llm_provided_name.split(f'-{parsed_stem}')[0]
        parsed_substituents = {}

        methyl_match = re.search(r'([\d,]+)-(\w*)methyl', substituent_part)
        if methyl_match:
            parsed_substituents["methyl"] = sorted([int(loc) for loc in methyl_match.group(1).split(',')])

        oxo_match = re.search(r'([\d,]+)-(\w*)oxo', substituent_part)
        if oxo_match:
            parsed_substituents["oxo"] = sorted([int(loc) for loc in oxo_match.group(1).split(',')])

        # 4. Check substituents
        if parsed_substituents != correct_properties["substituents"]:
            return f"Incorrect substituents. Expected {correct_properties['substituents']}, but parsed {parsed_substituents}."

        # 5. Check total carbon count
        calculated_carbons = parsed_length + len(parsed_substituents.get("methyl", []))
        if calculated_carbons != correct_properties["total_carbons"]:
            return f"Incorrect total carbon count. Expected {correct_properties['total_carbons']}, but calculated {calculated_carbons} from the name."

        # If all checks pass, the answer is correct.
        return "Correct"

    except Exception as e:
        return f"An error occurred during the checking process: {e}"

# Run the check and print the result
result = check_correctness_of_answer()
print(result)