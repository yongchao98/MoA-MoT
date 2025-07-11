def find_common_motifs():
    """
    Analyzes film motifs for Fritz Lang and William Friedkin to find a common element.
    """
    # Step 1 & 2: Store key imagery for each director, including film years.
    # We use lowercase keywords for easier matching.
    director_motifs = {
        "Fritz Lang": {
            "cyborg": "The first ever cyborg on screen (in 'Metropolis', 1927)",
            "bugs": "Bugs and spiders as symbols of madness/traps (in 'Dr. Mabuse the Gambler', 1922 and 'The Spiders', 1919)"
        },
        "William Friedkin": {
            "bugs": "Bugs as a central theme of paranoia and horror (in 'Bug', 2006, and locusts in 'The Exorcist', 1973)"
        }
    }

    # Step 3: Define the answer choices with keywords for matching.
    options = {
        "A": {"description": "Aboriginal masks", "keyword": "aboriginal masks"},
        "B": {"description": "Magic wands", "keyword": "magic wands"},
        "C": {"description": "The first ever cyborgs on screen", "keyword": "cyborg"},
        "D": {"description": "Bugs", "keyword": "bugs"}
    }

    print("Analyzing film motifs...\n")
    print(f"Fritz Lang's relevant motifs include: '{director_motifs['Fritz Lang']['cyborg']}' and '{director_motifs['Fritz Lang']['bugs']}'.")
    print(f"William Friedkin's relevant motifs include: '{director_motifs['William Friedkin']['bugs']}'.\n")
    
    final_answer = "E. None of the above"

    # Step 4 & 5: Iterate through options and check for commonality.
    for key, data in options.items():
        description = data["description"]
        keyword = data["keyword"]
        
        # Check if the keyword exists in both directors' motif dictionaries
        lang_has_motif = keyword in director_motifs["Fritz Lang"]
        friedkin_has_motif = keyword in director_motifs["William Friedkin"]
        
        is_common = "No"
        if lang_has_motif and friedkin_has_motif:
            is_common = "Yes"
            final_answer = f"{key}. {description}"

        print(f"Checking Option {key}: '{description}'")
        print(f"  - Present in Lang's work? {'Yes' if lang_has_motif else 'No'}")
        print(f"  - Present in Friedkin's work? {'Yes' if friedkin_has_motif else 'No'}")
        print(f"  - Is it a common motif? {is_common}\n")

    # Step 6: Print the conclusion
    print("---CONCLUSION---")
    print("The motif that appears in the works of both directors is:")
    print(final_answer)

find_common_motifs()