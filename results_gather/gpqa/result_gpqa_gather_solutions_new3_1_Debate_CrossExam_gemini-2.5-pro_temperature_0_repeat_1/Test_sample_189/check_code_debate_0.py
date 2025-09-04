import functools
import re

def check_nucleophile_order(candidate_answer_str: str) -> str:
    """
    Checks if the candidate's chosen option for nucleophile reactivity is correct.

    The function uses a rule-based system to determine the correct order of reactivity
    for the given nucleophiles in a polar protic solvent (aqueous solution).

    Nucleophiles:
    1. 4-methylcyclohexan-1-olate
    2. Hydroxide
    3. Propionate
    4. Methanol
    5. Ethanethiolate

    Args:
        candidate_answer_str: A string containing the candidate's answer, e.g., "<<<C>>>".

    Returns:
        A string indicating "Correct" or a detailed reason for the incorrectness.
    """

    # Define properties for each nucleophile
    # Steric hindrance scale: 1 (low), 2 (medium), 3 (high)
    nucleophiles_data = {
        1: {'name': '4-methylcyclohexan-1-olate', 'charge': -1, 'atom': 'O', 'resonance': False, 'steric': 3},
        2: {'name': 'Hydroxide', 'charge': -1, 'atom': 'O', 'resonance': False, 'steric': 1},
        3: {'name': 'Propionate', 'charge': -1, 'atom': 'O', 'resonance': True, 'steric': 2},
        4: {'name': 'Methanol', 'charge': 0, 'atom': 'O', 'resonance': False, 'steric': 1},
        5: {'name': 'Ethanethiolate', 'charge': -1, 'atom': 'S', 'resonance': False, 'steric': 2},
    }

    # Define periodic table period (higher number is lower in the table)
    atom_period = {'O': 2, 'S': 3}

    def compare_nucleophiles(nuc_id1, nuc_id2):
        """
        Compares two nucleophiles based on reactivity in a polar protic solvent.
        Returns -1 if nuc1 is more reactive, 1 if nuc2 is more reactive, 0 if equal.
        A more reactive nucleophile is "smaller" in sorting terms.
        """
        n1 = nucleophiles_data[nuc_id1]
        n2 = nucleophiles_data[nuc_id2]

        # Rule 1: Charge. Anions are much more reactive than neutral molecules.
        if n1['charge'] < 0 and n2['charge'] == 0:
            return -1  # n1 (anion) is more reactive
        if n2['charge'] < 0 and n1['charge'] == 0:
            return 1   # n2 (anion) is more reactive

        # Rule 2: Atom Identity (Polarizability in protic solvent).
        # Applies when both are anions.
        period1 = atom_period[n1['atom']]
        period2 = atom_period[n2['atom']]
        if period1 != period2:
            # Higher period number (lower in periodic table) is more nucleophilic in protic solvents.
            return -1 if period1 > period2 else 1

        # Rule 3: Resonance. Applies when both are anions with the same attacking atom.
        if n1['resonance'] and not n2['resonance']:
            return 1  # n2 (no resonance) is more reactive
        if n2['resonance'] and not n1['resonance']:
            return -1 # n1 (no resonance) is more reactive

        # Rule 4: Steric Hindrance. Applies when comparing non-resonant anions with the same atom.
        if n1['steric'] != n2['steric']:
            # Lower steric hindrance is more reactive.
            return -1 if n1['steric'] < n2['steric'] else 1

        return 0

    # Derive the correct order programmatically
    nucleophile_ids = list(nucleophiles_data.keys())
    # Use functools.cmp_to_key to sort using the custom comparison function
    correct_order = sorted(nucleophile_ids, key=functools.cmp_to_key(compare_nucleophiles))

    # Define the options from the question
    options = {
        'A': [2, 5, 1, 4, 3],
        'B': [5, 2, 3, 1, 4],
        'C': [5, 2, 1, 3, 4],
        'D': [2, 5, 3, 4, 3]
    }

    # Parse the candidate's answer
    match = re.search(r'<<<([A-D])>>>', candidate_answer_str)
    if not match:
        return "Incorrect format: The answer should be in the format <<<X>>> where X is A, B, C, or D."

    candidate_choice = match.group(1)

    if candidate_choice not in options:
        return f"Invalid option '{candidate_choice}' selected. Valid options are A, B, C, D."

    candidate_order = options[candidate_choice]

    # Check for invalid options like D, which has duplicates
    if len(set(candidate_order)) != 5:
        return f"The selected option '{candidate_choice}' with order {candidate_order} is invalid because it contains duplicate or missing numbers."

    # Compare the candidate's order with the correct order
    if candidate_order == correct_order:
        return "Correct"
    else:
        reason = f"The provided order {candidate_order} for option '{candidate_choice}' is incorrect.\n"
        reason += f"The correct order is {correct_order}, which corresponds to option C.\n\n"
        reason += "The reasoning is as follows:\n"
        reason += "1. **Ethanethiolate (5) is most reactive:** In a protic solvent (water), nucleophilicity increases down a group. Sulfur is larger and more polarizable than oxygen, and is less hindered by solvent hydrogen bonding, making it the strongest nucleophile.\n"
        reason += "2. **Methanol (4) is least reactive:** It is the only neutral molecule. Anions are always much stronger nucleophiles than their neutral counterparts.\n"
        reason += "3. **Propionate (3) is a weak nucleophile:** Its negative charge is delocalized by resonance across two oxygen atoms. This stabilization makes it much less reactive than anions with a localized charge.\n"
        reason += "4. **Hydroxide (2) vs. 4-methylcyclohexan-1-olate (1):** Hydroxide is more reactive. Although the alkoxide (1) is a slightly stronger base, it is significantly more sterically hindered by its bulky cyclohexyl ring. The small size of hydroxide (2) makes it a more effective nucleophile.\n"
        reason += f"Therefore, the correct ranking from most to least reactive is 5 > 2 > 1 > 3 > 4."
        return reason

# Example usage:
# print(check_nucleophile_order("<<<C>>>")) # Should return "Correct"
# print(check_nucleophile_order("<<<B>>>")) # Should return the reason for being incorrect
# print(check_nucleophile_order("<<<D>>>")) # Should return that the option is invalid