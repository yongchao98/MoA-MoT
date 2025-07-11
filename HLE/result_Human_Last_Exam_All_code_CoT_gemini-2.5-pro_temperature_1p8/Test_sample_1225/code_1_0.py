import textwrap

def find_equivalent_ballet_terms():
    """
    Analyzes pairs of ballet terms from the Royal Ballet School (RBS) and Vaganova Academy
    to find the equivalent pair.
    """
    ballet_knowledge = {
        'A': {
            "RBS": "Fifth position",
            "Vaganova": "Third position in arms",
            "is_equivalent": False,
            "explanation": "This pair is incorrect. 'Fifth position' refers to the feet, while 'third position in arms' refers to the arms. They are different parts of the body."
        },
        'B': {
            "RBS": "First arabesque",
            "Vaganova": "Third arabesque",
            "is_equivalent": True,
            "explanation": (
                "This is the correct choice. The equivalence is based on the fundamental principle of which arm is forward. "
                "In the English style (RBS), First Arabesque is defined by having the arm on the same side as the supporting leg extended forward. "
                "In the Vaganova system, both their First and Third arabesques follow this same rule (arm on the same side as the supporting leg is forward). "
                "While Vaganova's First and Third arabesques differ from each other in the height of the arms and torso alignment (épaulement), "
                "the Third Arabesque shares the same core defining principle as the English First Arabesque. The other options represent entirely different concepts or opposites."
            )
        },
        'C': {
            "RBS": "Assemblé",
            "Vaganova": "Brisé",
            "is_equivalent": False,
            "explanation": "This pair is incorrect. An assemblé is a jump to 'assemble' the feet, landing on two feet. A brisé is a 'broken', beaten jump that travels. They are two distinct steps."
        },
        'D': {
            "RBS": "Pirouette en dedan",
            "Vaganova": "Pirouette en dehor",
            "is_equivalent": False,
            "explanation": "This pair is incorrect. 'En dedan' means turning inward (towards the supporting leg), while 'en dehor' means turning outward (away from the supporting leg). They are opposites."
        },
        'E': {
            "RBS": "Temps levé",
            "Vaganova": "Sissonne",
            "is_equivalent": False,
            "explanation": "This pair is incorrect. A temps levé is a hop on one foot (taking off and landing on the same foot). A sissonne is a jump from two feet that lands on one foot. They have different mechanics."
        }
    }

    correct_option = None
    for option, details in ballet_knowledge.items():
        if details["is_equivalent"]:
            correct_option = option
            break
            
    if correct_option:
        details = ballet_knowledge[correct_option]
        print(f"The correct option is {correct_option}.")
        print(f"Term 1 (Royal Ballet School): {details['RBS']}")
        print(f"Term 2 (Vaganova Academy): {details['Vaganova']}")
        print("\nExplanation:")
        # Use textwrap for cleaner printing of long explanations
        wrapped_explanation = textwrap.fill(details['explanation'], width=80)
        print(wrapped_explanation)
    else:
        print("Could not determine the correct option.")

find_equivalent_ballet_terms()