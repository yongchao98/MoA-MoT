def solve_entailment_relation():
    """
    Determines the final projected natural logic operator for the given inference
    by breaking it down into logical steps and composing their relations.
    """
    # Define the 7 relations in MacCartney's framework
    RELATIONS = {
        '#': 'Forward Entailment',
        '>': 'Reverse Entailment',
        '≡': 'Equivalence',
        '^': 'Negation',
        '|': 'Alternation',
        'v': 'Cover',
        '/': 'Independence'
    }

    # Simplified composition table from MacCartney's framework
    # In the format: compose_table[relation1][relation2] = result_relation
    compose_table = {
        'Forward Entailment': {
            'Forward Entailment': 'Forward Entailment',
            # other compositions not needed for this problem
        },
        # other relations not needed for this problem
    }

    # 1. Define the logical steps
    premise = "Mark is singing a pop song by Taylor Swift"
    intermediate_premise = "Mark is singing a song by Taylor Swift"
    hypothesis = "Mark is not singing a song by Michael Jackson"

    print("Step-by-step logical transformation:\n")
    print(f"  P: \"{premise}\"")
    print(f"  P': \"{intermediate_premise}\"")
    print(f"  H: \"{hypothesis}\"\n")

    # 2. Determine the relation for the first step (P -> P')
    # A 'pop song' is a 'song', so P entails P'. This is a generalization.
    relation1 = 'Forward Entailment'
    print(f"Step 1: Relation between P and P'")
    print(f"  - The inference from P to P' is a generalization (a 'pop song' is a 'song').")
    print(f"  - This corresponds to the '{relation1}' operator.\n")


    # 3. Determine the relation for the second step (P' -> H)
    # Singing a song by one artist implies not singing one by another. P' entails H.
    relation2 = 'Forward Entailment'
    print(f"Step 2: Relation between P' and H")
    print(f"  - 'Singing a song by Taylor Swift' and 'singing a song by Michael Jackson' are mutually exclusive events.")
    print(f"  - Therefore, P' entails H.")
    print(f"  - This also corresponds to the '{relation2}' operator.\n")


    # 4. Compose the relations
    final_relation = compose_table[relation1][relation2]

    # 5. Print the final projected relation
    print("Final Projection:")
    print("The final relation is the composition of the relations from each step.")
    print("The equation is:")
    print(f"  {relation1} o {relation2} = {final_relation}")
    print(f"\nThe final projected natural logic operator is: {final_relation}")


solve_entailment_relation()
<<<Forward Entailment>>>