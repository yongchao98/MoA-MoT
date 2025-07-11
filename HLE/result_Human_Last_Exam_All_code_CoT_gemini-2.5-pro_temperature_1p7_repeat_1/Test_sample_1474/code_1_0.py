def find_common_cinematic_theme():
    """
    Analyzes the filmographies of Fritz Lang and William Friedkin to find a common
    theme from a given set of choices.
    """

    # Storing cinematic knowledge about the directors and themes.
    directors_data = {
        "Fritz Lang": {
            "A. Aboriginal masks": {
                "present": False,
                "evidence": "His films, while visually rich, do not have a known focus on Aboriginal masks."
            },
            "B. Magic wands": {
                "present": False,
                "evidence": "Themes of magic wands are not characteristic of Lang's expressionist, sci-fi, or noir films."
            },
            "C. The first ever cyborgs on screen": {
                "present": True,
                "evidence": "The iconic 'Maschinenmensch' (Machine-Person) in his 1927 film 'Metropolis' is widely considered a foundational cinematic cyborg/robot."
            },
            "D. Bugs": {
                "present": True,
                "evidence": "Insects appear as a motif in his work, notably in his early film serial 'The Spiders' ('Die Spinnen', 1919-20) and in hallucinatory sequences."
            }
        },
        "William Friedkin": {
            "A. Aboriginal masks": {
                "present": False,
                "evidence": "'The Exorcist' opens with the discovery of a Pazuzu amulet, an ancient Mesopotamian artifact, not an Aboriginal mask."
            },
            "B. Magic wands": {
                "present": False,
                "evidence": "Friedkin's realistic and often gritty style does not incorporate fantasy elements like magic wands."
            },
            "C. The first ever cyborgs on screen": {
                "present": False,
                "evidence": "Friedkin's work does not feature cyborgs."
            },
            "D. Bugs": {
                "present": True,
                "evidence": "Bugs are a significant theme. The title of his 2006 film is 'Bug', which centers on a paranoia about insect infestation. 'The Exorcist' also features the unsettling 'spider-walk' scene."
            }
        }
    }

    print("Investigating shared themes between Fritz Lang and William Friedkin...")
    print("-" * 60)

    final_answer_letter = "E"
    common_theme_description = "None of the above"

    choices = ["A. Aboriginal masks", "B. Magic wands", "C. The first ever cyborgs on screen", "D. Bugs"]

    for choice in choices:
        lang_analysis = directors_data["Fritz Lang"][choice]
        friedkin_analysis = directors_data["William Friedkin"][choice]

        print(f"Analyzing Choice: {choice}")
        print(f"  - For Fritz Lang: {'Present' if lang_analysis['present'] else 'Absent'}. {lang_analysis['evidence']}")
        print(f"  - For William Friedkin: {'Present' if friedkin_analysis['present'] else 'Absent'}. {friedkin_analysis['evidence']}")

        if lang_analysis["present"] and friedkin_analysis["present"]:
            print(f"  - Verdict: This theme is present in both directors' oeuvres.\n")
            final_answer_letter = choice[0]
            common_theme_description = choice
        else:
            print(f"  - Verdict: This theme is NOT common to both.\n")

    print("-" * 60)
    print("Conclusion:")
    print(f"The common theme found in the work of both Fritz Lang and William Friedkin is: {common_theme_description}")

# Run the analysis
find_common_cinematic_theme()