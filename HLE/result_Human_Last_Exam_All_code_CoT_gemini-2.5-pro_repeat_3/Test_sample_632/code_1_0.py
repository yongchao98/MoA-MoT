def find_historical_feature():
    """
    This function simulates a database lookup to find a historical
    feature of a famous Roman square.
    """
    # A small knowledge base of Roman piazzas and their historical features.
    piazza_knowledge_base = {
        "Piazza della Rotonda": {
            "current_main_attraction": "The Pantheon",
            "feature_removed_around_1950": "A bustling and historic fish market (known as La Pescheria)"
        },
        "Piazza Navona": {
            "current_main_attraction": "Fountain of the Four Rivers",
            "feature_removed_around_1950": "The main food market of Rome (moved in 1869)"
        }
    }

    target_piazza = "Piazza della Rotonda"

    # Retrieve the information from the knowledge base
    if target_piazza in piazza_knowledge_base:
        feature_info = piazza_knowledge_base[target_piazza]
        removed_feature = feature_info.get("feature_removed_around_1950")

        if removed_feature:
            print(f"The unique feature of the '{target_piazza}' that was removed was:")
            print(removed_feature)
        else:
            print(f"No removed feature is listed for '{target_piazza}'.")
    else:
        print(f"Information about '{target_piazza}' not found.")

find_historical_feature()

<<<A fish market>>>