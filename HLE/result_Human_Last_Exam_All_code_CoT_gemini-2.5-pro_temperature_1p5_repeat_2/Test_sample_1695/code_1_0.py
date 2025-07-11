def solve_sentence_puzzle():
    """
    Solves the sentence permutation puzzle by defining subjects, predicates,
    and a set of rules to evaluate the truthfulness of each combination.
    """

    # 1. Define the components of the sentences
    subjects = ["Satyendranath Bose", "The conversation", "Sushi", "Our Tuna"]
    predicates = [
        ("baptized", "Bosons"),
        ("turned", "playful"),
        ("contains", "fish"),
        ("is", "a Bengalese")
    ]

    # 2. Establish a "World Model" with properties for each subject
    subject_properties = {
        "Satyendranath Bose": {"type": "person", "origin": "Bengalese"},
        "The conversation": {"type": "abstract"},
        "Sushi": {"type": "food"},
        "Our Tuna": {"type": "animal", "identity": "Bengalese"} # e.g., Bengal cat
    }

    # 3. Define a function to evaluate if a sentence is true
    def is_true(subject, verb, obj):
        props = subject_properties[subject]
        subj_type = props["type"]

        # Predicate 1: "baptized Bosons"
        if verb == "baptized" and obj == "Bosons":
            # Only a person with historical significance in this context can do this.
            return subject == "Satyendranath Bose"

        # Predicate 2: "turned playful"
        if verb == "turned" and obj == "playful":
            # Things with behavior can become playful.
            return subj_type in ["person", "animal", "abstract"]

        # Predicate 3: "contains fish"
        if verb == "contains" and obj == "fish":
            # Foods can be made of ingredients; animals can eat food.
            return subj_type in ["food", "animal"]

        # Predicate 4: "is a Bengalese"
        if verb == "is" and obj == "a Bengalese":
            # Identity check for origin or breed.
            if subj_type == "person":
                return props.get("origin") == "Bengalese"
            if subj_type == "animal":
                return props.get("identity") == "Bengalese"
            return False

        return False

    # 4. Iterate, count, and prepare the equation
    counts_per_subject = []
    for s in subjects:
        count = 0
        for v, o in predicates:
            if is_true(s, v, o):
                count += 1
        counts_per_subject.append(count)

    # 5. Print the final calculation and result
    equation = " + ".join(map(str, counts_per_subject))
    total = sum(counts_per_subject)
    print(f"{equation} = {total}")

solve_sentence_puzzle()
<<<8>>>