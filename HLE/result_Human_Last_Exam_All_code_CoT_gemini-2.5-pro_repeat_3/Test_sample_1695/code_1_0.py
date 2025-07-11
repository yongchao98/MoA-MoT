def find_true_and_meaningful_sentences():
    """
    This function finds all "true and meaningful" sentences that can be formed
    by combining subjects, verbs, and objects from the four initial sentences.
    """
    subjects = ['Satyendranath Bose', 'The conversation', 'Sushi', 'Our Tuna']
    objects = ['Bosons', 'playful', 'fish', 'Bengalese']
    
    # Sentence structures (verb phrases)
    structures = {
        'baptized': 1,
        'turned': 2,
        'contains': 3,
        'is a': 4
    }

    true_sentences = []
    counts_per_structure = {1: 0, 2: 0, 3: 0, 4: 0}

    print("--- Finding True and Meaningful Sentences ---")

    for s in subjects:
        for o in objects:
            # Check Structure 1: baptized
            if s == 'Satyendranath Bose' and o == 'Bosons':
                sentence = f"{s} baptized {o}."
                if sentence not in true_sentences:
                    true_sentences.append(sentence)
                    counts_per_structure[1] += 1
            
            # Check Structure 2: turned
            if o == 'playful' and (s == 'The conversation' or s == 'Our Tuna'):
                sentence = f"{s} turned {o}."
                if sentence not in true_sentences:
                    true_sentences.append(sentence)
                    counts_per_structure[2] += 1
            
            # Check Structure 3: contains
            if s == 'Sushi' and o == 'fish':
                sentence = f"{s} contains {o}."
                if sentence not in true_sentences:
                    true_sentences.append(sentence)
                    counts_per_structure[3] += 1

            # Check Structure 4: is a
            if o == 'Bengalese' and (s == 'Satyendranath Bose' or s == 'Our Tuna'):
                sentence = f"{s} is a {o}."
                if sentence not in true_sentences:
                    true_sentences.append(sentence)
                    counts_per_structure[4] += 1

    # Print the discovered sentences
    for sentence in sorted(true_sentences):
        print(sentence)
    
    print("\n--- Calculation ---")
    c1 = counts_per_structure[1]
    c2 = counts_per_structure[2]
    c3 = counts_per_structure[3]
    c4 = counts_per_structure[4]
    total = c1 + c2 + c3 + c4

    # Print the equation as requested
    print(f"Number of true sentences per structure (baptized, turned, contains, is a):")
    print(f"{c1} + {c2} + {c3} + {c4} = {total}")
    
    # Final answer in the required format
    print(f"\nTotal number of true and meaningful sentences is {total}.")


if __name__ == '__main__':
    find_true_and_meaningful_sentences()
    # The final answer is 6
    print("\n<<<6>>>")
