def compare_covalency(overlap_comparison):
    """
    Determines the relative covalency based on orbital overlap.

    Args:
        overlap_comparison (str): A string describing the relative orbital overlap
                                  ('greater', 'weaker', 'equal').

    Returns:
        str: The corresponding relative covalency ('stronger', 'weaker', 'equal').
    """
    if overlap_comparison == "greater":
        return "stronger"
    elif overlap_comparison == "weaker":
        return "weaker"
    else:
        return "equal"

# Information from the problem statement:
# The 4f orbital overlap in CeF6(2-) is greater than in CeCl6(2-).
overlap_in_CeF6_vs_CeCl6 = "greater"

# Covalency is directly proportional to orbital overlap.
# A greater overlap results in stronger covalency.
covalency_of_CeF6 = compare_covalency(overlap_in_CeF6_vs_CeCl6)

# Print the reasoning and the result.
# There are no numerical equations in this conceptual problem.
print("Step 1: The strength of covalency is directly related to the extent of orbital overlap.")
print("Step 2: The problem states that the orbital overlap in CeF6(2-) is 'greater' than in CeCl6(2-).")
print(f"Step 3: Therefore, because it has greater orbital overlap, CeF6(2-) will exhibit '{covalency_of_CeF6}' covalency compared to CeCl6(2-).")

<<<stronger>>>