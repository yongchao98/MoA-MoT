def find_french_equivalent():
    """
    Finds and prints the French equivalent of a given Gaulish name based on a predefined mapping.
    """
    # A dictionary mapping known Gaulish/Latinized names to their French equivalents.
    # The name Odobeccus comes from a Germanic root *aud-, meaning 'wealth',
    # which gave rise to the Latin name Odo and the French name Eudes.
    gaulish_to_french_map = {
        "Odobeccus": "Eudes"
    }

    # The name we are looking up
    gaulish_name = "Odobeccus"

    # Retrieve the corresponding French name
    if gaulish_name in gaulish_to_french_map:
        french_name = gaulish_to_french_map[gaulish_name]
        print(f"The French equivalent of the Gaulish name '{gaulish_name}' is '{french_name}'.")
    else:
        print(f"Could not find a French equivalent for '{gaulish_name}'.")

# Execute the function to get the answer
find_french_equivalent()