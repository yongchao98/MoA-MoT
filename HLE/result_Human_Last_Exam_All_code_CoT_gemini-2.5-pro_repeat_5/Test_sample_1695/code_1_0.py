def solve_sentence_puzzle():
    """
    This script solves the sentence puzzle by systematically generating and
    evaluating all possible sentences based on a defined set of semantic rules.
    """

    # Step 1: Define entities and their properties (our knowledge base)
    subjects_data = {
        "Satyendranath Bose": {'type': 'person', 'name': 'Satyendranath Bose'},
        "The conversation": {'type': 'abstract', 'name': 'The conversation'},
        "Sushi": {'type': 'food', 'name': 'Sushi'},
        "Our Tuna": {'type': 'animal', 'animal_type': 'cat', 'name': 'Our Tuna'}
    }

    objects_data = {
        "Bosons": {'type': 'physics_concept', 'name': 'Bosons'},
        "playful": {'type': 'attribute', 'name': 'playful'},
        "fish": {'type': 'food', 'name': 'fish'},
        "a Bengalese": {'type': 'animal_descriptor', 'animal_type': 'cat', 'name': 'a Bengalese'}
    }

    subjects = list(subjects_data.keys())
    objects = list(objects_data.keys())

    # Step 2: Define the rules for each verb as a function
    def is_true_baptized(subj_props, obj_props):
        """Rule for the verb 'baptized'."""
        # Only true for the specific scientific naming context.
        return subj_props['name'] == 'Satyendranath Bose' and obj_props['name'] == 'Bosons'

    def is_true_turned(subj_props, obj_props):
        """Rule for the verb 'turned'."""
        # A subject can 'turn' into a state or attribute.
        if obj_props['type'] == 'attribute' and obj_props['name'] == 'playful':
            # Persons, animals, and abstract concepts can become playful.
            return subj_props['type'] in ['person', 'animal', 'abstract']
        return False

    def is_true_contains(subj_props, obj_props):
        """Rule for the verb 'contains'."""
        # True in the context of food ingredients.
        return subj_props['name'] == 'Sushi' and obj_props['name'] == 'fish'

    def is_true_is(subj_props, obj_props):
        """Rule for the verb 'is'/'is a'."""
        # Case 1: Identity/type description ('Our Tuna is a Bengalese')
        if obj_props['type'] == 'animal_descriptor':
            return subj_props.get('animal_type') == obj_props.get('animal_type')

        # Case 2: Attribute description ('X is playful')
        if obj_props['type'] == 'attribute' and obj_props['name'] == 'playful':
            return subj_props['type'] in ['person', 'animal', 'abstract']

        # Case 3: Composition description ('Sushi is fish')
        if subj_props['name'] == 'Sushi' and obj_props['name'] == 'fish':
            return True
            
        return False

    verbs_and_rules = [
        ("baptized", is_true_baptized),
        ("turned", is_true_turned),
        ("contains", is_true_contains),
        ("is", is_true_is)
    ]

    # Step 3 & 4: Iterate, evaluate all sentences, and print results
    true_sentences = []
    total_count = 0

    print("The following sentences are considered true and meaningful:")

    sum_parts = []
    for verb, rule_function in verbs_and_rules:
        group_count = 0
        sentences_in_group = []
        for subj_name in subjects:
            for obj_name in objects:
                subj_props = subjects_data[subj_name]
                obj_props = objects_data[obj_name]
                
                if rule_function(subj_props, obj_props):
                    total_count += 1
                    group_count += 1
                    # The original sentence for 'is' is 'is a'. Let's reconstruct it.
                    verb_display = "is a" if verb == "is" and obj_props['type'] == 'animal_descriptor' else verb
                    sentence = f"- {subj_name} {verb_display} {obj_name}."
                    sentences_in_group.append(sentence)
        
        if group_count > 0:
            sum_parts.append(str(group_count))
            for s in sorted(sentences_in_group):
                print(s)

    print("\nCalculating the total number of true and meaningful sentences:")
    equation = " + ".join(sum_parts)
    print(f"Total = {equation} = {total_count}")

solve_sentence_puzzle()