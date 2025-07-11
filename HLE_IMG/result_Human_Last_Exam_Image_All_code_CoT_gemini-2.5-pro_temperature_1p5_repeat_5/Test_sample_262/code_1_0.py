def analyze_complex_stability():
    """
    Analyzes the stability of four iridium complexes based on their structure.

    The stability of these emitters in LECs is largely determined by the
    presence of fluorine atoms ortho to the C-Ir bond. This position is
    known to be a weak point leading to molecular degradation and shorter
    device lifetimes.
    """

    # We represent the presence of an ortho-fluorine as a boolean.
    # True = ortho-fluorine is present (unstable, shorter lifetime)
    # False = no ortho-fluorine (stable, longer lifetime)
    complexes_properties = {
        1: {'has_ortho_fluorine': False, 'description': "Fluorines are meta and para to the C-Ir bond."},
        2: {'has_ortho_fluorine': False, 'description': "Ligands have no fluorine or non-ortho fluorines."},
        3: {'has_ortho_fluorine': True,  'description': "Contains a fluorine atom ortho to the C-Ir bond."},
        4: {'has_ortho_fluorine': True,  'description': "Contains a fluorine atom ortho to the C-Ir bond."}
    }

    shorter_lifetime_complexes = []
    print("Analysis of Complexes:")
    for i, properties in complexes_properties.items():
        stability = "shorter" if properties['has_ortho_fluorine'] else "longer"
        print(f"Complex {i}: {properties['description']} -> Expected to have a {stability} lifetime.")
        if properties['has_ortho_fluorine']:
            shorter_lifetime_complexes.append(i)

    # Sort the list for consistent output, e.g., [3, 4] not [4, 3]
    shorter_lifetime_complexes.sort()
    
    print("\nConclusion:")
    print(f"The complexes expected to show shorter lifetimes are those with fluorine atoms ortho to the C-Ir bond.")
    # The formatting requires printing each number in the final list
    # The required answer format is [3, 4], which corresponds to option J.
    # To be clear, let's print the list.
    print(f"This corresponds to complexes: {shorter_lifetime_complexes}")
    
    # Mapping the list to the answer choices
    answer_map = {
        '[1]': 'A', '[2]': 'B', '[3]': 'C', '[4]': 'D',
        '[1, 2]': 'E', '[1, 3]': 'F', '[1, 4]': 'G',
        '[2, 3]': 'H', '[2, 4]': 'I', '[3, 4]': 'J',
        '[1, 2, 3]': 'K', '[1, 2, 4]': 'L', '[1, 3, 4]': 'M',
        '[2, 3, 4]': 'N', '[1, 2, 3, 4]': 'O'
    }
    
    answer_key = str(shorter_lifetime_complexes)
    final_answer_letter = answer_map.get(answer_key, "Unknown")
    
    print(f"The correct answer choice is {final_answer_letter}.")


analyze_complex_stability()
# The final answer is the combination of complexes 3 and 4.
# Based on the options, this is choice J.
# The code should output J as the answer
# but the instruction said "use 'print' function for the output when relevant"
# and "directly return the answer with the format <<<answer content>>>"
# which implies I should determine the answer and then print it at the end.
# So I should print J.

final_answer = 'J'
print(f'<<<J>>>')
