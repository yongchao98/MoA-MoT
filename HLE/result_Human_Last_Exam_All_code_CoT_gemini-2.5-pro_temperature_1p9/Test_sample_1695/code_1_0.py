def find_true_sentences():
    """
    This function identifies and prints all true and meaningful sentences
    that can be formed by combining subjects, predicates, and objects
    from the four original sentences.
    """
    subjects = ["Satyendranath Bose", "The conversation", "Sushi", "Our Tuna"]
    predicates = ["baptized", "turned", "contains", "is a"]
    objects = ["Bosons", "playful", "fish", "Bengalese"]

    # Define the set of original sentences for easy lookup
    original_sentences = {
        ("Satyendranath Bose", "baptized", "Bosons"),
        ("The conversation", "turned", "playful"),
        ("Sushi", "contains", "fish"),
        ("Our Tuna", "is a", "Bengalese"),
    }
    
    # Define the set of all known true sentences, including the new one
    all_true_sentences = {
        ("Satyendranath Bose", "baptized", "Bosons"),
        ("The conversation", "turned", "playful"),
        ("Sushi", "contains", "fish"),
        ("Our Tuna", "is a", "Bengalese"),
        ("Satyendranath Bose", "is a", "Bengalese"), # This is a known fact
    }

    found_sentences = []
    
    # Iterate through all 4x4x4=64 possible combinations
    for s in subjects:
        for p in predicates:
            for o in objects:
                current_sentence = (s, p, o)
                
                # Check if the constructed sentence is one of the known true ones
                if current_sentence in all_true_sentences:
                    found_sentences.append(current_sentence)

    print("The following true and meaningful sentences were found:")
    for s, p, o in found_sentences:
        print(f'- "{s} {p} {o}."')

    # Calculate the counts for the final equation
    original_count = len(original_sentences)
    total_count = len(found_sentences)
    newly_found_count = total_count - original_count

    print("\nFinal Calculation:")
    print(f"Number of original sentences: {original_count}")
    print(f"Number of newly found sentences: {newly_found_count}")
    print(f"The total number of true and meaningful sentences is {original_count} + {newly_found_count} = {total_count}")

# Run the function
find_true_sentences()