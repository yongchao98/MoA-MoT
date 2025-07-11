def identify_opera():
    """
    This script identifies an opera by matching clues from its score
    against a small knowledge base.
    """
    # --- Clues extracted from the provided sheet music ---
    # The language for instruments and directions is German.
    # e.g., "gr. Fl." (große Flöten), "engl. Horn", "(Vorhang)", "(mit Dämpfer)"
    language_clue = "German"

    # The score includes a specific stage direction at the beginning of the act.
    stage_direction_clue = "(Vorhang)"  # German for "(Curtain)"

    # The prompt specifies this is the opening of the second act.
    context_clue = "Opening of Act 2"
    
    # The lush, large orchestra points to a specific musical era.
    style_clue = "Late Romantic / Early 20th Century"

    # --- Simulated Knowledge Base ---
    opera_database = {
        "Der Rosenkavalier": {
            "composer": "Richard Strauss",
            "language": "German",
            "style": "Late Romantic / Early 20th Century",
            "act_2_opening_notes": "The curtain rises ('Vorhang') immediately on a scene of great activity."
        },
        "Die Meistersinger von Nürnberg": {
            "composer": "Richard Wagner",
            "language": "German",
            "style": "Romantic",
            "act_2_opening_notes": "Opens with an evening scene on a street in Nuremberg."
        },
        "The Marriage of Figaro": {
            "composer": "Wolfgang Amadeus Mozart",
            "language": "Italian",
            "style": "Classical",
            "act_2_opening_notes": "Opens in the Countess's boudoir."
        }
    }

    # --- Matching Logic ---
    identified_opera = "Unknown"
    for opera, data in opera_database.items():
        # Check if the clues match the data for the current opera
        is_language_match = (data["language"] == language_clue)
        is_style_match = (style_clue in data["style"])
        is_stage_direction_match = (stage_direction_clue in data["act_2_opening_notes"])

        if is_language_match and is_style_match and is_stage_direction_match:
            identified_opera = opera
            break

    print(f"Based on the clues, the opera is identified as: {identified_opera}")

identify_opera()