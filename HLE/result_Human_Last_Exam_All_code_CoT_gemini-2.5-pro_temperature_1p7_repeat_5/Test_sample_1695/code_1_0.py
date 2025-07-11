def solve_sentence_puzzle():
    """
    This script finds all true and meaningful sentences by combining the subjects,
    verbs, and objects from the four provided sentences. It then calculates the
    total count and displays it as a sum.
    """
    
    # Step 1: Define the components of the original sentences
    subjects = ["Satyendranath Bose", "The conversation", "Sushi", "Our Tuna"]
    predicates = ["baptized", "turned", "contains", "is a"]
    objects = ["Bosons", "playful", "fish", "Bengalese"]

    # Step 2: Establish a set to hold the valid sentences found
    true_sentences = set()

    # Step 3: Apply logical and factual rules to find all true sentences
    # Rule for "baptized": Requires a person as a subject.
    # "Satyendranath Bose baptized Bosons" is historically true.
    true_sentences.add("Satyendranath Bose baptized Bosons.")

    # Rule for "turned": Subject changes state to an adjective ('playful').
    # "The conversation turned playful" is true (original).
    true_sentences.add("The conversation turned playful.")
    # "Our Tuna turned playful" is true (assuming Tuna is a cat/pet).
    true_sentences.add("Our Tuna turned playful.")

    # Rule for "contains": An entity having a component.
    # "Sushi contains fish" is true (original).
    true_sentences.add("Sushi contains fish.")
    # "Our Tuna contains fish" is true (a cat that has eaten fish).
    true_sentences.add("Our Tuna contains fish.")

    # Rule for "is a": Subject belongs to a class (noun). The phrase "is a"
    # cannot be followed by an adjective like "playful".
    # "Our Tuna is a Bengalese" is true (original).
    true_sentences.add("Our Tuna is a Bengalese.")
    # "Satyendranath Bose is a Bengalese" is true (he was from Bengal).
    true_sentences.add("Satyendranath Bose is a Bengalese.")

    # Step 4: Count the number of true sentences for each predicate
    counts_per_predicate = {p: 0 for p in predicates}
    for sentence in true_sentences:
        # This loop correctly attributes each found sentence to its verb.
        words = sentence.split()
        verb = words[1]
        if verb == "is": # Handles the "is a" case
            verb = "is a"
        counts_per_predicate[verb] +=1
    
    counts = list(counts_per_predicate.values())
    total_count = sum(counts)

    # Step 5: Print the result as a final equation
    print(f"{counts[0]} + {counts[1]} + {counts[2]} + {counts[3]} = {total_count}")

solve_sentence_puzzle()
<<<7>>>