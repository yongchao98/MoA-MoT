import textwrap

def solve_avalanche():
    """
    Analyzes the Word Avalanche problem and prints the correct answer with justification.
    """
    description = "My software tells the birds when and where to relieve themselves."
    keyword = "computer"
    choices = {
        'A': "Computers comp few urns",
        'B': "Computer: \"Come poo, tern!\"",
        'C': "Computer: \"Come, pee you turkey!\"",
        'D': "Comb pewter computer",
        'E': "Computer needs bird relieves"
    }

    print("# Task Description")
    print(textwrap.fill(f"Find a true avalanche using the word '{keyword}' that fits the description: '{description}'\n", 80))

    print("# Analysis")
    print(textwrap.fill("A true avalanche requires a phonetic repetition of syllables and must fit the given context. The key word is 'computer', which has the syllables 'com-pu-ter'.\n", 80))

    # Detailed analysis of the best choice
    print("## Evaluating the Best Choice (B)\n")
    choice_b_text = choices['B']
    print(f"Choice B is: {choice_b_text}\n")

    print("# 1. Phonetic Repetition Analysis")
    print(textwrap.fill("The phrase 'Come poo, tern!' is a strong phonetic repetition of 'com-pu-ter'.", 80))
    print("- 'Come' sounds like 'com'.")
    print("- 'poo' sounds like 'pu'.")
    print("- 'tern' sounds like 'ter'.\n")


    print("# 2. Contextual Fit Analysis")
    print(textwrap.fill(f"The phrase depicts a 'Computer' telling a 'tern' (a type of sea bird) to 'poo' (relieve itself). This perfectly matches the description: '{description}'.\n", 80))

    print("# Conclusion")
    print(textwrap.fill("Choice B is the only option that is both a strong phonetic pun and a perfect match for the context described.", 80))
    print("-" * 30)
    print(f"The correct answer is B: {choice_b_text}")
    print("<<<B>>>")


solve_avalanche()