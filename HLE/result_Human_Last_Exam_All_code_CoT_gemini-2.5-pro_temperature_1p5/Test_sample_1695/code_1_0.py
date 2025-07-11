def find_true_sentences():
    """
    This function forms all possible sentences by combining subjects and predicates,
    evaluates their truthfulness and meaning, and prints the results.
    """
    subjects = [
        "Satyendranath Bose",
        "The conversation",
        "Sushi",
        "Our Tuna"
    ]

    predicates = [
        "baptized Bosons.",
        "turned playful.",
        "contains fish.",
        "is a Bengalese."
    ]

    # We pre-determine which combinations are "true and meaningful" based on world knowledge.
    # The original sentences are (0,0), (1,1), (2,2), (3,3).
    # Other plausible combinations:
    # - "Our Tuna turned playful." -> subject 3, predicate 1 -> True, cats are playful.
    # - "Our Tuna contains fish." -> subject 3, predicate 2 -> True, cat food is often fish-based.
    true_combinations = {(0, 0), (1, 1), (2, 2), (3, 3), (3, 1), (3, 2)}

    print("The true and meaningful sentences are:")
    
    true_sentences_list = []
    original_count = 0
    new_count = 0

    # Iterate through all 16 possible combinations
    for i, subject in enumerate(subjects):
        for j, predicate in enumerate(predicates):
            if (i, j) in true_combinations:
                sentence = f"{subject} {predicate}"
                true_sentences_list.append(sentence)
                if i == j:
                    original_count += 1
                else:
                    new_count += 1
    
    for sentence in true_sentences_list:
        print(f"- {sentence}")
        
    print("\n---")
    print("Final Count Calculation:")
    print(f"Original true sentences = {original_count}")
    print(f"Newly formed true sentences = {new_count}")
    print(f"Total number of true and meaningful sentences = {original_count} + {new_count} = {len(true_sentences_list)}")

find_true_sentences()