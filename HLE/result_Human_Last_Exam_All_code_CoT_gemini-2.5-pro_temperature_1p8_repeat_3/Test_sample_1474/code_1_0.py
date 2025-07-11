def find_common_movie_themes():
    """
    This script analyzes themes in the films of Fritz Lang and William Friedkin
    to find a common element from a given list of options.
    """
    
    # Step 1: Define notable themes/imagery for each director based on film history.
    # For Lang, "Bugs" references the hallucinatory insect swarms in "Dr. Mabuse the Gambler" (1922).
    # For Friedkin, "Bugs" is the central topic of his psychological horror film "Bug" (2006).
    lang_themes = {
        "The first ever cyborgs on screen",  # From "Metropolis"
        "Bugs"                               # From "Dr. Mabuse the Gambler"
    }

    friedkin_themes = {
        "Demonic possession", # From "The Exorcist"
        "Bugs"                # From "Bug"
    }

    # Step 2: Define the answer choices.
    options = {
        "A": "Aboriginal masks",
        "B": "Magic wands",
        "C": "The first ever cyborgs on screen",
        "D": "Bugs",
        "E": "None of the above"
    }

    print("Analyzing themes for Fritz Lang and William Friedkin...\n")
    print("Checking each option:")

    correct_option_key = "E"  # Default to "E" if no commonality is found

    # Step 3: Iterate through options A-D to check for commonality.
    for key in ["A", "B", "C", "D"]:
        theme = options[key]
        in_lang = theme in lang_themes
        in_friedkin = theme in friedkin_themes
        
        print(f"- Checking Option {key} ('{theme}'):")
        print(f"  - Appears in Lang's work? {'Yes' if in_lang else 'No'}")
        print(f"  - Appears in Friedkin's work? {'Yes' if in_friedkin else 'No'}")
        
        if in_lang and in_friedkin:
            correct_option_key = key
            print("  -> Found a common element!")
            break
        else:
            print("  -> Not a common element.")
        print()

    print("\n--- CONCLUSION ---")

    # Step 4: Print the "final equation" showing the intersection of themes.
    common_themes_found = lang_themes.intersection(friedkin_themes)
    
    # Outputting each element as requested.
    lang_list = list(lang_themes)
    friedkin_list = list(friedkin_themes)
    common_list = list(common_themes_found)

    print("The final equation represents the intersection of the two directors' themes:")
    print(f"{lang_list[0]}, {lang_list[1]}  (Lang)  ∩  {friedkin_list[0]}, {friedkin_list[1]}  (Friedkin)  =  {common_list[0]}")

    print(f"\nThe correct option is '{correct_option_key}' because '{options[correct_option_key]}' appears in the works of both directors.")

find_common_movie_themes()