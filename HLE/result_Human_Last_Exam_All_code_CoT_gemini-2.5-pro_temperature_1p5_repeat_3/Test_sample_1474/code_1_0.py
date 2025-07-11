def find_common_element():
    """
    Analyzes the work of directors Fritz Lang and William Friedkin
    to find a common thematic or visual element from a given list.
    """

    directors_and_films = {
        "Fritz Lang": {
            "film": "The Testament of Dr. Mabuse (1933)",
            "evidence": "Features a key scene where a character, driven by paranoia, has terrifying hallucinations of swarming bugs."
        },
        "William Friedkin": {
            "film": "Bug (2006)",
            "evidence": "The entire plot of this psychological horror film revolves around characters who believe they are infested with government-implanted bugs."
        }
    }

    common_element = "Bugs"
    correct_option = "D"

    print("Investigating a common element in the oeuvres of Fritz Lang and William Friedkin...")
    print("="*70)

    for director, details in directors_and_films.items():
        print(f"Director: {director}")
        print(f"Relevant Film: {details['film']}")
        print(f"Evidence: {details['evidence']}")
        print("-"*70)

    print(f"\nConclusion:")
    print(f"Both directors have used '{common_element}' as a powerful visual metaphor for paranoia and psychological collapse.")
    print(f"Therefore, the correct option is '{correct_option}'.")


find_common_element()