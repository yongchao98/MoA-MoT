import functools

def check_answer():
    """
    Checks the correctness of the provided LLM answer by programmatically ranking the nucleophiles
    based on established chemical principles for a polar protic solvent.
    """

    # 1. Define the nucleophiles with their key chemical properties for ranking.
    # Steric hindrance is simplified to 'low' vs 'high' for comparison.
    nucleophiles = [
        {'id': 1, 'name': '4-methylcyclohexan-1-olate', 'charge': 'anion', 'atom': 'O', 'resonance': False, 'sterics': 'high'},
        {'id': 2, 'name': 'Hydroxide', 'charge': 'anion', 'atom': 'O', 'resonance': False, 'sterics': 'low'},
        {'id': 3, 'name': 'Propionate', 'charge': 'anion', 'atom': 'O', 'resonance': True, 'sterics': 'low'},
        {'id': 4, 'name': 'Methanol', 'charge': 'neutral', 'atom': 'O', 'resonance': False, 'sterics': 'low'},
        {'id': 5, 'name': 'Ethanethiolate', 'charge': 'anion', 'atom': 'S', 'resonance': False, 'sterics': 'low'},
    ]

    # 2. Create a comparison function that applies the rules of nucleophilicity in a protic solvent.
    # The function returns -1 if nuc1 > nuc2, 1 if nuc2 > nuc1, and 0 if they are equal.
    def compare_nucleophiles(nuc1, nuc2):
        # Rule 1: Charge (Anion > Neutral)
        if nuc1['charge'] == 'anion' and nuc2['charge'] == 'neutral':
            return -1  # nuc1 is stronger
        if nuc1['charge'] == 'neutral' and nuc2['charge'] == 'anion':
            return 1   # nuc2 is stronger

        # If charges are the same (both anions), apply next rules.
        # Rule 2: Attacking Atom (S > O in protic solvent due to polarizability/solvation)
        if nuc1['atom'] == 'S' and nuc2['atom'] == 'O':
            return -1  # nuc1 is stronger
        if nuc1['atom'] == 'O' and nuc2['atom'] == 'S':
            return 1   # nuc2 is stronger

        # If atoms are the same (both O), apply next rules.
        # Rule 3: Resonance (No Resonance > Resonance)
        if not nuc1['resonance'] and nuc2['resonance']:
            return -1  # nuc1 is stronger
        if nuc1['resonance'] and not nuc2['resonance']:
            return 1   # nuc2 is stronger

        # If resonance is the same (both have no resonance), apply next rule.
        # Rule 4: Steric Hindrance (Low > High)
        steric_map = {'low': 0, 'high': 1}
        if steric_map[nuc1['sterics']] < steric_map[nuc2['sterics']]:
            return -1  # nuc1 is stronger (less hindered)
        if steric_map[nuc1['sterics']] > steric_map[nuc2['sterics']]:
            return 1   # nuc2 is stronger (less hindered)

        return 0 # Should not be reached in this specific problem set

    # 3. Sort the list of nucleophiles to determine the chemically correct order.
    sorted_nucleophiles = sorted(nucleophiles, key=functools.cmp_to_key(compare_nucleophiles))
    correct_order = [n['id'] for n in sorted_nucleophiles]

    # 4. Extract the order and final choice from the provided answer.
    # The provided answer's reasoning gives the order 5, 2, 1, 3, 4.
    # The final choice is <<<B>>>.
    llm_reasoning_order = [5, 2, 1, 3, 4]
    llm_final_choice = 'B'

    # 5. Define the options based on the candidate answers to verify the final choice.
    # Option B is consistently shown as 5, 2, 1, 3, 4.
    options = {
        'A': [5, 2, 3, 1, 4],
        'B': [5, 2, 1, 3, 4],
    }

    # 6. Check the correctness of the provided answer.
    # Check 1: Is the order derived in the reasoning chemically correct?
    if llm_reasoning_order != correct_order:
        return (f"Incorrect. The reasoning in the answer leads to a chemically incorrect order of nucleophilicity.\n"
                f"Expected order based on principles: {correct_order}\n"
                f"Answer's derived order: {llm_reasoning_order}")

    # Check 2: Does the final letter choice (e.g., <<<B>>>) correspond to the correct order?
    if options.get(llm_final_choice) != correct_order:
        return (f"Incorrect. The final choice '{llm_final_choice}' does not match the correct chemical order.\n"
                f"The correct order is {correct_order}, which corresponds to option B in the provided context.\n"
                f"The answer incorrectly selected option {llm_final_choice}.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result.
print(check_answer())