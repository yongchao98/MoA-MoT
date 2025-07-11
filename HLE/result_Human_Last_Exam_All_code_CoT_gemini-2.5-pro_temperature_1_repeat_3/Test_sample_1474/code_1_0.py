def solve_movie_trivia():
    """
    Analyzes common themes in the films of Fritz Lang and William Friedkin
    to find the correct answer from a list of choices.
    """
    
    # This dictionary stores the analysis for each director and theme.
    # The format is: (Does the theme appear?, "Justification")
    analysis = {
        "Fritz Lang": {
            "A. Aboriginal masks": (False, "No prominent use of Aboriginal masks in Lang's major films."),
            "B. Magic wands": (False, "Lang's films focus on crime, dystopia, and fate, not fantasy wands."),
            "C. The first ever cyborgs on screen": (True, "The 'Maschinenmensch' in 'Metropolis' (1927) is a famous proto-cyborg/automaton."),
            "D. Bugs": (True, "Lang directed the two-part silent film 'Die Spinnen' ('The Spiders') in 1919-1920, where the antagonists are a spider-themed criminal syndicate.")
        },
        "William Friedkin": {
            "A. Aboriginal masks": (False, "'The Exorcist' features a Mesopotamian (not Aboriginal) demon artifact."),
            "B. Magic wands": (False, "Friedkin's style is grounded in realism or visceral horror, not fantasy."),
            "C. The first ever cyborgs on screen": (False, "Friedkin's filmography does not feature this theme, and he did not direct the first film with cyborgs."),
            "D. Bugs": (True, "Friedkin directed the psychological horror film 'Bug' (2006). Furthermore, 'The Exorcist' features the infamous 'spider walk' scene and the demon Pazuzu's association with locusts.")
        }
    }

    choices = [
        "A. Aboriginal masks",
        "B. Magic wands",
        "C. The first ever cyborgs on screen",
        "D. Bugs"
    ]

    correct_answer = "E. None of the above"
    final_choice_letter = "E"
    explanation_lang = ""
    explanation_friedkin = ""

    # Iterate through choices to find one that is true for both directors
    for choice in choices:
        lang_has_theme = analysis["Fritz Lang"][choice][0]
        friedkin_has_theme = analysis["William Friedkin"][choice][0]

        if lang_has_theme and friedkin_has_theme:
            correct_answer = choice
            final_choice_letter = choice.split('.')[0]
            explanation_lang = analysis["Fritz Lang"][choice][1]
            explanation_friedkin = analysis["William Friedkin"][choice][1]
            break

    print("Step-by-step analysis of common themes:")
    print("="*40)
    print(f"Checking Choice C: 'The first ever cyborgs on screen'")
    print(f"- Fritz Lang: YES. {analysis['Fritz Lang']['C. The first ever cyborgs on screen'][1]}")
    print(f"- William Friedkin: NO. {analysis['William Friedkin']['C. The first ever cyborgs on screen'][1]}")
    print("Conclusion: Choice C is not common to both.")
    print("-" * 40)
    print(f"Checking Choice D: 'Bugs'")
    print(f"- Fritz Lang: YES. {explanation_lang}")
    print(f"- William Friedkin: YES. {explanation_friedkin}")
    print("Conclusion: Choice D is a theme/image found in the works of both directors.")
    print("="*40)
    print(f"\nThe correct shared theme is: {correct_answer}")
    
    print(f"\n<<<{final_choice_letter}>>>")

solve_movie_trivia()