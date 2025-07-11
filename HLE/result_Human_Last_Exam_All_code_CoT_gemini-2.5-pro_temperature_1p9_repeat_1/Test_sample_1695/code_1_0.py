def solve_sentence_puzzle():
    """
    This function finds and counts the number of true and meaningful sentences
    that can be formed by swapping subjects and objects from a given set of
    four sentences.
    """
    subjects = ["Satyendranath Bose", "The conversation", "Sushi", "Our Tuna"]
    verbs = ["baptized", "turned", "contains", "is"]
    objects = ["Bosons", "playful", "fish", "a Bengalese"]

    # This list contains tuples of (subject, verb, object) for all sentences
    # determined to be true and meaningful based on real-world knowledge.
    true_combinations = [
        # Verb: baptized (action of naming, typically by a person)
        # Only Satyendranath Bose (the physicist) "baptizing" (naming) Bosons makes sense.
        ("Satyendranath Bose", "baptized", "Bosons"),

        # Verb: turned (to become a certain state/quality)
        # An animate object (person, cat) or an abstract concept (conversation) can turn playful.
        ("Satyendranath Bose", "turned", "playful"),
        ("The conversation", "turned", "playful"),
        ("Our Tuna", "turned", "playful"),

        # Verb: contains (to have something inside)
        # Sushi is made with fish, and a cat can contain fish it has eaten.
        ("Sushi", "contains", "fish"),
        ("Our Tuna", "contains", "fish"),

        # Verb: is (describes a state of being or identity)
        # A person, a cat, or a conversation can be described as playful.
        # Sushi can be equated with fish as its primary ingredient.
        # The cat Tuna is identified as a Bengalese.
        ("Satyendranath Bose", "is", "playful"),
        ("The conversation", "is", "playful"),
        ("Sushi", "is", "fish"),
        ("Our Tuna", "is", "playful"),
        ("Our Tuna", "is", "a Bengalese"),
    ]

    # Calculate counts for each verb structure
    counts = {verb: 0 for verb in verbs}
    for sub, verb, obj in true_combinations:
        if verb in counts:
            counts[verb] += 1

    total_count = len(true_combinations)
    
    # Create the equation string as requested
    count_values = [counts[verb] for verb in verbs]
    equation_str = " + ".join(map(str, count_values))

    print("The final count is the sum of true sentences for each verb structure:")
    print(f"{equation_str} = {total_count}")

solve_sentence_puzzle()