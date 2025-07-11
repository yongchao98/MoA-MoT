import textwrap

def analyze_liability():
    """
    This function analyzes the legal scenario and determines the correct allocation of liability.
    """

    # Analysis of Mark's Incident
    mark_analysis = """
    Mark was an employee of Evergreen Grass Care Ltd. acting within the scope of his employment.
    His actions (getting distracted, playing with a dog while on a mower) were negligent.
    This negligence directly caused damage to the pool.
    Liability: Mark is personally liable for his negligence. Evergreen Grass Care Ltd. is vicariously liable for the actions of its employee.
    Conclusion for Mark: Evergreen and Mark are jointly and severally liable.
    """

    # Analysis of Lincoln's Incident
    lincoln_analysis = """
    Lincoln was also an employee acting within the scope of his employment.
    His action (using a blower on rocks near a car) was negligent, causing foreseeable damage.
    The extent of the damage (minor scratches) affects the amount of compensation owed but does not eliminate liability.
    Liability: Lincoln is personally liable for his negligence. Evergreen Grass Care Ltd. is vicariously liable.
    Conclusion for Lincoln: Evergreen and Lincoln are jointly and severally liable.
    """

    # Final Conclusion based on the answer choices
    final_conclusion = "The correct statement is that Evergreen and Mark are jointly and severally liable for the damage from Mark's actions, AND Evergreen and Lincoln are jointly and severally liable for the damage from Lincoln's actions."

    # Map conclusion to the answer choices
    answer_choice = "E"
    answer_text = "Evergreen Grass Care Ltd. and Mark are jointly and severally liable for the damage that resulted from Mark's actions. Evergreen Grass Care Ltd. and Lincoln are jointly and severally liable for the damage that resulted from Lincoln's actions."

    print("Analyzing the legal scenario to determine liability...")
    print("-" * 50)
    print("Incident 1: Mark and the Pool")
    print(textwrap.dedent(mark_analysis))
    print("-" * 50)
    print("Incident 2: Lincoln and the Ferrari")
    print(textwrap.dedent(lincoln_analysis))
    print("-" * 50)
    print("Final Conclusion:")
    print(final_conclusion)
    print(f"\nThis corresponds to Answer Choice: {answer_choice}")
    print(f'"{answer_text}"')


if __name__ == '__main__':
    analyze_liability()