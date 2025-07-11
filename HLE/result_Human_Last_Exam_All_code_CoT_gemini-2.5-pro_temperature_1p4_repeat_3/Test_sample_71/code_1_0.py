import re

def solve_stress_puzzle():
    """
    This script determines the stressed syllable in Old Russian phrases
    by applying a derived set of linguistic rules.
    """

    # The list of phrases we need to analyze
    phrases_to_analyze = [
        "i ne znali",
        "i povelo že",
        "ne vymyla že",
        "ponesla",
        "vyvela že",
        "i unesli"
    ]

    # This dictionary stores the results based on the rule-based analysis.
    # The logic for each is explained in the comments.
    results = {
        # 'ne' + non-feminine -> stress on 'ne'. Syllables: i(1)-ne(2)-zna(3)-li(4). Result: 2
        "i ne znali": 2,
        # 'po-' prefix -> stress on root 've'. Syllables: i(1)-po(2)-ve(3)-lo(4)-že(5). Result: 3
        "i povelo že": 3,
        # 'vy-' prefix -> stress on 'vy'. Syllables: ne(1)-vy(2)-my(3)-la(4)-že(5). Result: 2
        "ne vymyla že": 2,
        # 'po-' prefix -> stress on root 'nes'. Syllables: po(1)-ne(2)-sla(3). Result: 2
        "ponesla": 2,
        # 'vy-' prefix -> stress on 'vy'. Syllables: vy(1)-ve(2)-la(3)-že(4). Result: 1
        "vyvela že": 1,
        # 'u-' prefix -> stress on root 'nes'. Syllables: i(1)-u(2)-ne(3)-sli(4). Result: 3
        "i unesli": 3,
    }

    final_answer_string = ""
    print("Analysis of each phrase:")
    print("="*30)

    for phrase in phrases_to_analyze:
        stress_syllable_number = results[phrase]
        final_answer_string += str(stress_syllable_number)
        
        # To show the "equation", we print each component number.
        print(f"For '{phrase}':")
        print(f"The stressed syllable is number {stress_syllable_number}.")
        print("-" * 30)

    print(f"\nThe combined six-digit result is: {final_answer_string}")

solve_stress_puzzle()
<<<232213>>>