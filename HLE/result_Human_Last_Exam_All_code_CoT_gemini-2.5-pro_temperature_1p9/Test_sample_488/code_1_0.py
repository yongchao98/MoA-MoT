import textwrap

def analyze_metaphor():
    """
    Analyzes the conceptual metaphor in the phrase "my love for humanity"
    based on the work of Lakoff and Johnson.
    """

    phrase = "my love for humanity"
    choices = {
        'A': 'Conceptual metaphor',
        'B': 'Physical metaphor',
        'C': 'Orientational metaphor',
        'D': 'Structural metaphor',
        'E': 'Intentional metaphor'
    }

    print(f"Analyzing the conceptual metaphor in the phrase: '{phrase}'\n")

    print("Step 1: Identify the underlying metaphorical mapping.")
    explanation_1 = ("The phrase treats 'love,' which is an abstract emotion, as if it were a concrete "
                     "substance or object that can be possessed. The use of 'my' indicates possession, "
                     "similar to 'my book' or 'my coffee'. This establishes a metaphorical mapping: "
                     "LOVE IS A POSSESSABLE SUBSTANCE/OBJECT.")
    print(textwrap.fill(explanation_1, width=80))
    print("-" * 20)

    print("Step 2: Evaluate the given answer choices against this mapping.")
    explanation_c = ("An Orientational Metaphor gives concepts a spatial orientation (like up-down, in-out). "
                     "For example, 'HAPPY IS UP' ('I'm feeling up.'). This does not apply to '{phrase}'.")
    explanation_d = ("A Structural Metaphor occurs when one concept is metaphorically structured in terms of another. "
                     "By conceptualizing love as an object or substance, we give it a structure—it can be "
                     "measured, contained, and, as in this case, possessed. This fits the mapping perfectly.")
    
    print(f"Choice C ({choices['C']}):")
    print(textwrap.fill(explanation_c.format(phrase=phrase), width=80))
    print()
    print(f"Choice D ({choices['D']}):")
    print(textwrap.fill(explanation_d, width=80))
    print("-" * 20)
    
    print("Step 3: Reach a conclusion.")
    conclusion = ("The most specific term for this is an 'ontological metaphor,' which deals with viewing "
                  "abstract concepts as physical entities. Since that is not an option, we choose the best "
                  "available category. The metaphor provides a structure (the structure of an object) to "
                  "the concept of love. Therefore, it is best classified as a Structural Metaphor.")
    print(textwrap.fill(conclusion, width=80))

    final_answer = 'D'
    print("\nFinal Answer Choice:", final_answer)


if __name__ == "__main__":
    analyze_metaphor()

<<<D>>>