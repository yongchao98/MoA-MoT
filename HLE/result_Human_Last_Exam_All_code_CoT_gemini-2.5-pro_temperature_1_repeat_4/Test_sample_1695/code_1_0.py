import collections

def solve_sentence_puzzle():
    """
    Solves the sentence swapping puzzle by identifying and counting
    all true and meaningful combinations.
    """
    subjects = ["Satyendranath Bose", "The conversation", "Sushi", "Our Tuna"]
    verbs = ["baptized", "turned", "contains", "is a"]
    objects = ["Bosons", "playful", "fish", "Bengalese"]

    # Based on analysis, the following combinations of (subject_idx, verb_idx, object_idx)
    # are deemed true and meaningful.
    # Interpretation: "Our Tuna" is a Bengal cat. "Satyendranath Bose" was from Bengal.
    true_sentences_indices = [
        (0, 0, 0),  # Satyendranath Bose baptized Bosons. (Metaphorically true)
        (1, 1, 1),  # The conversation turned playful. (Given)
        (2, 2, 2),  # Sushi contains fish. (Generally true)
        (3, 3, 3),  # Our Tuna is a Bengalese. (Given, establishes Tuna as a Bengal cat)
        (0, 3, 3),  # Satyendranath Bose is a Bengalese. (Factually true)
        (3, 2, 2)   # Our Tuna contains fish. (True for a cat's diet)
    ]

    print("The true and meaningful sentences are:")
    
    # Sort for consistent output, originals first, then new ones
    # Original indices are (0,0,0), (1,1,1), (2,2,2), (3,3,3)
    # We can sort by tuple to group them logically
    true_sentences_indices.sort()
    
    verb_counts = collections.defaultdict(int)
    count = 0

    for i, (sub_idx, verb_idx, obj_idx) in enumerate(true_sentences_indices):
        count += 1
        sentence = f"{subjects[sub_idx]} {verbs[verb_idx]} {objects[obj_idx]}."
        print(f"{i+1}. {sentence}")
        verb_counts[verbs[verb_idx]] += 1

    print("\nThe number of true sentences for each verb category are:")
    calc_parts = []
    for verb in verbs:
        num = verb_counts[verb]
        print(f"- {verb}: {num}")
        calc_parts.append(str(num))

    equation = " + ".join(calc_parts)
    total = sum(verb_counts.values())

    print(f"\nFinal calculation: {equation} = {total}")
    print(f"Total number of true and meaningful sentences: {total}")

solve_sentence_puzzle()
<<<6>>>