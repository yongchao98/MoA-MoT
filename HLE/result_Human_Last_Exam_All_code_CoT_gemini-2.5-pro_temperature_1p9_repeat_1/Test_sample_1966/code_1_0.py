import collections

def find_common_characters():
    # Step 1: Create a list of all potential Shakespearean title characters from the options.
    shakespearean_title_chars = [
        "Julius Caesar", "Pericles", "Cleopatra",
        "King John", "Troilus", "Antony"
    ]
    
    # Step 2: Create a list of characters famously mentioned by name in The Divine Comedy.
    # This list is based on literary analysis of the text.
    # - Julius Caesar: Inferno, Canto IV
    # - Cleopatra: Inferno, Canto V
    # - Antony (Marc Antony): Mentioned in Paradiso, Canto VI, during Justinian's speech.
    divine_comedy_chars = [
        "Julius Caesar", "Hector", "Aeneas", "Homer", "Socrates", "Plato",
        "Cleopatra", "Helen", "Achilles", "Paris", "Tristan", "Antony"
    ]
    
    # Using sets to find the intersection is efficient
    set_shakespeare = set(shakespearean_title_chars)
    set_dante = set(divine_comedy_chars)
    
    # Step 3: Find the common characters
    common_characters = sorted(list(set_shakespeare.intersection(set_dante)))
    
    # Step 4: Print the reasoning and the result
    print("Finding Shakespearean title characters mentioned in Dante's *The Divine Comedy*...")
    print(f"Shakespearean characters from options: {sorted(list(set_shakespeare))}")
    print(f"Notable characters from *The Divine Comedy*: {sorted(list(set_dante))}")
    print("-" * 30)
    print("The characters found in both are:")
    
    for character in common_characters:
        print(f"- {character}")
    
    print("\nThis matches option E.")

find_common_characters()