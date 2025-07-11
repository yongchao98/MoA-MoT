import textwrap

def solve_micromalthidae_question():
    """
    This function explains the reasoning to answer a question about the feeding habits
    of an adult male Micromalthidae beetle.
    """
    question = "Suppose an adult male Micromalthidae is found in the colony. Upon its death, what will be the only thing that this individual will have fed on?"
    
    choices = {
        'A': 'Its mother',
        'B': 'Fungus',
        'C': 'Decaying wood',
        'D': 'Cellulose',
        'E': 'Nothing'
    }

    print("Analyzing the biological question:")
    print("-" * 30)
    print(textwrap.fill(question, width=80))
    print("-" * 30)

    # Step-by-step reasoning
    print("\nReasoning:\n")
    
    explanation = [
        "1. The organism in question is an adult male beetle from the family Micromalthidae.",
        "2. While the larvae of Micromalthidae feed on decaying wood, the adult stage is very different.",
        "3. Scientific literature confirms that adult male Micromalthidae are a rare, short-lived, terminal stage.",
        "4. Most importantly, they possess vestigial (non-functional) mouthparts. They are physically incapable of eating.",
        "5. Therefore, during its entire life as an adult, the male beetle does not feed on anything.",
        "6. Although the beetle fed on wood as a larva, the question focuses on the adult individual. The most distinct feeding characteristic of the adult male is the complete absence of feeding.",
        "7. Thus, upon its death, the adult male will have consumed nothing during its adult life stage."
    ]
    
    for step in explanation:
        print(textwrap.fill(step, width=80))

    # Stating the final answer
    correct_choice = 'E'
    print("\n" + "-" * 30)
    print(f"Conclusion: Based on the biology of the adult male, the correct answer is E.")
    print(f"The adult male beetle feeds on '{choices[correct_choice]}'.")
    print("-" * 30)


solve_micromalthidae_question()