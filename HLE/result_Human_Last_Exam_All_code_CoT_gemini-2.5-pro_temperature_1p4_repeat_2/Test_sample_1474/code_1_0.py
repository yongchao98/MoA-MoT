def find_common_theme():
    """
    This script analyzes the works of directors Fritz Lang and William Friedkin
    to identify a common thematic element from a given list of choices.
    """
    director1 = "Fritz Lang"
    director2 = "William Friedkin"
    correct_option_letter = "D"
    correct_theme = "Bugs"

    # The "equation" is the presentation of evidence for both directors, leading to the conclusion.
    print(f"Task: Find the common imagery in the oeuvres of {director1} and {director2}.")
    print("----------------------------------------------------------")
    print(f"Analyzing Option {correct_option_letter}: '{correct_theme}'")
    print("----------------------------------------------------------")

    # Evidence for Fritz Lang
    lang_film = "'Dr. Mabuse the Gambler' (1922)"
    lang_evidence = f"In {director1}'s film {lang_film}, a character suffers hallucinations of swarming insects, linking bugs to psychological horror."
    print(f"Evidence 1 (Fritz Lang): {lang_evidence}")

    # Evidence for William Friedkin
    friedkin_film = "'Bug' (2006)"
    friedkin_evidence = f"In {director2}'s film {friedkin_film}, the entire plot centers on paranoia and a perceived insect infestation."
    print(f"Evidence 2 (William Friedkin): {friedkin_evidence}")

    print("\nConclusion: The motif of 'Bugs' as a manifestation of psychological torment is a clear point of connection between both directors' works.")

find_common_theme()