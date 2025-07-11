import collections

def solve_counseling_case():
    """
    Analyzes the counseling options for a mother whose adolescent son is vaping.
    """
    # Step 1: Define the statements and their validity based on medical principles.
    # True = Valid advice, False = Invalid advice
    statements = {
        'I': {
            "text": "Vaping is a better option for her son than cigarettes. It is ok for him to continue vaping, as long as he is trying to cut down.",
            "is_valid": False,
            "reason": "This normalizes nicotine use in an adolescent. The goal should be complete cessation, not harm reduction by switching products."
        },
        'II': {
            "text": "It would be good for her son to start using nicotine patches, gum, or lozenges instead of vaping.",
            "is_valid": True,
            "reason": "Nicotine Replacement Therapy (NRT) is a valid, evidence-based first-line treatment for nicotine addiction to aid in cessation."
        },
        'III': {
            "text": "Vaping’s risks and benefits are well understood in adults, but not in children, so her son should not vape at all.",
            "is_valid": True,
            "reason": "This is a crucial educational point. It correctly highlights the unknown long-term risks for adolescents and advises complete cessation, directly addressing the mother's perspective."
        },
        'IV': {
            "text": "Vaping’s risks and benefits remain poorly understood in children, but it has been shown vaping has clear benefits over cigarettes in children.",
            "is_valid": False,
            "reason": "This language is dangerous. Framing any nicotine product as having 'clear benefits' for a child is misleading and inappropriate."
        },
        'V': {
            "text": "Consider initiating bupropion and varenicline depending on her son’s needs.",
            "is_valid": True,
            "reason": "These are valid, evidence-based prescription medication options for adolescent nicotine cessation."
        }
    }

    # Step 2: Define the multiple-choice answers.
    AnswerChoice = collections.namedtuple('AnswerChoice', ['label', 'components'])
    answer_choices = [
        AnswerChoice('A', ['I']), AnswerChoice('B', ['II']), AnswerChoice('C', ['III']),
        AnswerChoice('D', ['IV']), AnswerChoice('E', ['V']), AnswerChoice('F', ['I', 'II']),
        AnswerChoice('G', ['I', 'III']), AnswerChoice('H', ['I', 'IV']), AnswerChoice('I', ['I', 'V']),
        AnswerChoice('J', ['II', 'III']), AnswerChoice('K', ['II', 'IV']), AnswerChoice('L', ['II', 'V']),
        AnswerChoice('M', ['III', 'IV']), AnswerChoice('N', ['III', 'V']), AnswerChoice('O', ['IV', 'V']),
        AnswerChoice('P', ['I', 'II', 'III']), AnswerChoice('Q', ['II', 'III', 'IV']), AnswerChoice('R', ['I', 'III', 'IV']),
        AnswerChoice('S', ['I', 'II', 'IV']), AnswerChoice('T', ['III', 'IV', 'V']), AnswerChoice('U', ['I', 'IV', 'V']),
        AnswerChoice('V', ['II', 'IV', 'V'])
    ]

    # Step 3: Filter for choices that ONLY contain valid statements.
    valid_choices = []
    print("Evaluating all possible answers...\n")
    for choice in answer_choices:
        is_choice_valid = all(statements[comp]['is_valid'] for comp in choice.components)
        if is_choice_valid:
            valid_choices.append(choice)
            print(f"Answer choice {choice.label} ({', '.join(choice.components)}) contains only valid statements.")

    print("\n--- Applying further reasoning to find the BEST answer ---\n")

    # Step 4: Determine the best choice among the valid ones based on counseling principles.
    # The mother's positive experience with vaping must be addressed. Statement III does this best.
    best_choice = None
    reasoning = ""
    for choice in valid_choices:
        if 'III' in choice.components:
            # Between J(II, III) and N(III, V), J is superior because NRT (II) is a more common
            # first-line recommendation than prescription medications (V).
            # The combination of the core rationale (III) and the most common first-line treatment (II)
            # is the most robust initial counseling strategy.
            if 'II' in choice.components:
                 best_choice = choice
                 reasoning = (
                     "The best counseling strategy must include Statement III to educate the mother\n"
                     "and correct her misconception that vaping is safe for her son.\n"
                     "Among the valid options including III, the combination of III and II is optimal.\n"
                     "It pairs the essential rationale (Statement III) with a concrete, accessible,\n"
                     "and standard first-line treatment option (Statement II: NRT)."
                 )
                 break # Found the best choice

    print(f"Final Reasoning: {reasoning}\n")
    print(f"The best combination of options is: {', '.join(best_choice.components)}")
    print(f"This corresponds to answer choice: {best_choice.label}")
    print("\nThe chosen options are:")
    for comp in best_choice.components:
        print(f"Statement {comp}: {statements[comp]['text']}")

    final_answer = f"<<<{best_choice.label}>>>"
    print(f"\nFinal Answer: {final_answer}")

solve_counseling_case()