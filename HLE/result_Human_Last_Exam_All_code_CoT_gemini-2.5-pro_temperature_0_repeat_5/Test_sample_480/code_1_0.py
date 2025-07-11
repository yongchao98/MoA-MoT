def solve():
    """
    Calculates the final projected natural logic operator for the given inference.
    """
    # The 7 relations: Equivalence, Forward Entailment, Reverse Entailment,
    # Negation, Alternation, Cover, Independence
    relations = ["≡", "<", ">", "!", "|", "^", "⌄"]
    rel_map = {name: i for i, name in enumerate(relations)}
    
    # MacCartney's composition table B o A
    # Rows are B (context), columns are A (edit)
    composition_table = [
        # ≡   <   >   !   |   ^   ⌄
        [ 0,  1,  2,  3,  4,  5,  6 ], # ≡
        [ 1,  1,  6,  4,  4,  6,  6 ], # <
        [ 2,  6,  2,  5,  6,  5,  6 ], # >
        [ 3,  2,  1,  0,  5,  4,  6 ], # !
        [ 4,  4,  6,  5,  5,  6,  6 ], # |
        [ 5,  6,  5,  4,  6,  4,  6 ], # ^
        [ 6,  6,  6,  6,  6,  6,  6 ]  # ⌄
    ]

    def compose(b_idx, a_idx):
        return composition_table[b_idx][a_idx]

    # Step 1: Define atomic edits and their contexts
    edits = [
        {"name": "is singing -> is not singing", "atomic": "!", "context": "≡"},
        {"name": "pop -> ε (pop song -> song)", "atomic": "<", "context": "!"},
        {"name": "Taylor Swift -> Michael Jackson", "atomic": "|", "context": "!"}
    ]
    
    print("Premise: 'Mark is singing a pop song by Taylor Swift'")
    print("Hypothesis: 'Mark is not singing a song by Michael Jackson'")
    print("-" * 20)
    print("Analyzing edits from left to right:")
    
    # Step 2: Calculate projected relations
    projected_relations = []
    for edit in edits:
        context_idx = rel_map[edit["context"]]
        atomic_idx = rel_map[edit["atomic"]]
        proj_idx = compose(context_idx, atomic_idx)
        projected_relations.append(proj_idx)
        print(f"Edit: '{edit['name']}'")
        print(f"  - Atomic Relation: {edit['atomic']}")
        print(f"  - Context Relation: {edit['context']}")
        print(f"  - Projected Relation ({edit['context']} ∘ {edit['atomic']}): {relations[proj_idx]}")

    # Step 3: Compose projected relations sequentially
    print("-" * 20)
    print("Composing projected relations sequentially:")
    
    final_rel_idx = projected_relations[0]
    print(f"Start with first projected relation: {relations[final_rel_idx]}")
    
    # Composition: (Proj1 o Proj2) o Proj3
    # First composition
    res_tmp_idx = compose(projected_relations[0], projected_relations[1])
    print(f"({relations[projected_relations[0]]} ∘ {relations[projected_relations[1]]}) = {relations[res_tmp_idx]}")
    
    # Final composition
    final_rel_idx = compose(res_tmp_idx, projected_relations[2])
    print(f"({relations[res_tmp_idx]} ∘ {relations[projected_relations[2]]}) = {relations[final_rel_idx]}")

    final_relation_name = {
        "≡": "Equivalence",
        "<": "Forward Entailment",
        ">": "Reverse Entailment",
        "!": "Negation",
        "|": "Alternation",
        "^": "Cover",
        "⌄": "Independence"
    }[relations[final_rel_idx]]
    
    print("-" * 20)
    print(f"The final projected operator is {relations[final_rel_idx]}, which is named '{final_relation_name}'.")

solve()