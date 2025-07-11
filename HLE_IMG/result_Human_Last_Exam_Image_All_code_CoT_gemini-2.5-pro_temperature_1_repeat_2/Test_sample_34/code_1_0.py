import textwrap

def analyze_kinship_diagram():
    """
    This function interprets the Lévi-Strauss kinship diagram, derives its structural pattern,
    and evaluates the given anthropological systems to find the correct match.
    """
    print("Step 1: Decoding the symbols in the Lévi-Strauss diagram.")
    print("-" * 60)
    print("Δ = Male\no = Female\n= = Marriage (Husband/Wife)\nHorizontal line = Siblings (Brother/Sister)\nDescending line = Descent (Parent/Child)\n+ = Familiar/Close relationship\n- = Reserved/Distant relationship")
    print("-" * 60, "\n")

    print("Step 2: Deriving the 'kinship equation' from the specific diagram.")
    print("-" * 60)
    
    # This section outputs the 'numbers' (i.e., the +/- signs) for the 'equation'.
    print("The four key relationships and their attitudes are:")
    print("1. Husband / Wife: The '=' has a '-' sign. The relationship is Reserved (-).")
    print("2. Brother / Sister: The horizontal line has a '+' sign. The relationship is Familiar (+).")
    print("3. Father / Son: The descending line from the father has a '+' sign. The relationship is Familiar (+).")
    print("4. Mother's Brother / Sister's Son: The line to the son from the uncle has a '-' sign. The relationship is Reserved (-).")
    print("-" * 60, "\n")

    print("Step 3: Interpreting the overall system type.")
    print("-" * 60)
    explanation = """
    The pattern is: Fa/So (+), MB/ZS (-). This indicates a system where the paternal line is dominant and the bond is strong, while the maternal uncle's role is distant. This is the hallmark of a PATRILINEAL system.
    """
    print(textwrap.dedent(explanation).strip())
    print("-" * 60, "\n")

    print("Step 4: Evaluating the answer choices.")
    print("-" * 60)
    evaluation = """
    - A, B, E: Incorrect. These include matrilineal societies (Trobriand, Siuoi), which have a different structural emphasis, typically the inverse of the diagram for male relationships (strong maternal uncle role).
    - D: Incorrect. Cherkess and Tonga are classically described with the opposite pattern: a reserved Father/Son (-) and a familiar Mother's Brother/Son (+) relationship.
    - C: Correct. This option contains two patrilineal systems. The Lake Kubutu system is known to match the diagram's pattern. By this match and the elimination of other options, this is the most plausible answer.
    """
    print(textwrap.dedent(evaluation).strip())
    print("-" * 60)

# Run the analysis
analyze_kinship_diagram()

<<<C>>>