import textwrap

def analyze_liability():
    """
    This function analyzes the provided legal scenario and prints the correct answer choice
    along with the justification.
    """

    # The scenario involves two separate incidents of negligence by employees of Evergreen Grass Care Ltd.
    # We must apply the principle of vicarious liability, where an employer is responsible
    # for the negligent acts of its employees performed within the scope of their employment.

    # Incident 1: Mark and the pool
    # Mark is directly liable for his negligence.
    # Evergreen is vicariously liable for Mark's actions.
    liability_mark = "Evergreen Grass Care Ltd. and Mark are jointly and severally liable for the damage that resulted from Mark's actions."

    # Incident 2: Lincoln and the Ferrari
    # Lincoln is directly liable for his negligence. The "minimal" nature of the damage affects the
    # amount of compensation, not the existence of liability itself.
    # Evergreen is vicariously liable for Lincoln's actions.
    liability_lincoln = "Evergreen Grass Care Ltd. and Lincoln are jointly and severally liable for the damage that resulted from Lincoln's actions."

    # Combining the two conclusions gives us the correct answer.
    correct_answer_choice = "E"
    full_answer_text = f"{liability_mark}  {liability_lincoln}"

    print(f"The correct answer is: {correct_answer_choice}")
    print("\n" + "="*20 + " JUSTIFICATION " + "="*20)
    print("\n".join(textwrap.wrap(full_answer_text, width=70)))
    print("\n" + "="*55)


analyze_liability()