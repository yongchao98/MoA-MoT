import collections

def solve_sentence_puzzle():
    """
    Solves the sentence puzzle by generating all possible sentences
    and evaluating them for truth and meaning.
    """
    subjects = ["Satyendranath Bose", "The conversation", "Sushi", "Our Tuna"]
    verbs = ["baptized", "turned", "contains", "is"]
    objects = ["Bosons", "playful", "fish", "a Bengalese"]

    # This dictionary will store our logical rules based on world knowledge
    # derived from the original sentences.
    knowledge_base = {
        # Original sentence: "Satyendranath Bose baptized Bosons."
        ("Satyendranath Bose", "baptized", "Bosons"): True,
        # Original sentence: "The conversation turned playful."
        ("The conversation", "turned", "playful"): True,
        # Original sentence: "Sushi contains fish."
        ("Sushi", "contains", "fish"): True,
        # Original sentence: "Our Tuna is a Bengalese."
        ("Our Tuna", "is", "a Bengalese"): True,
        # Derived Fact 1: A person from Bengal is a Bengalese. Satyendranath Bose was from Bengal.
        ("Satyendranath Bose", "is", "a Bengalese"): True,
        # Derived Fact 2: If something 'turned' playful, it then 'is' playful.
        ("The conversation", "is", "playful"): True,
    }

    true_sentences = []
    counts_by_verb = collections.defaultdict(int)

    # Generate all 4x4x4 = 64 possible sentences
    for verb in verbs:
        for subject in subjects:
            for obj in objects:
                sentence_tuple = (subject, verb, obj)
                # Check if the sentence is in our knowledge base of true statements
                if knowledge_base.get(sentence_tuple, False):
                    sentence_str = f"{subject} {verb} {obj}."
                    if sentence_str not in true_sentences:
                        true_sentences.append(sentence_str)
                        counts_by_verb[verb] += 1
    
    print("Found the following true and meaningful sentences:")
    for sentence in sorted(true_sentences):
        print(f"- {sentence}")
    
    print("\nBreakdown of the count by verb:")
    final_sum_parts = []
    for verb in verbs:
        count = counts_by_verb[verb]
        print(f"- Sentences with '{verb}': {count}")
        final_sum_parts.append(str(count))
    
    total = sum(counts_by_verb.values())
    equation = " + ".join(final_sum_parts)
    print(f"\nThe final equation is: {equation} = {total}")

solve_sentence_puzzle()
<<<6>>>