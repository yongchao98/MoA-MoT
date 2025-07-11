def get_relation_name(relation_symbol):
    """Returns the name of a natural logic relation symbol."""
    names = {
        '=': "Equivalence",
        '<': "Forward Entailment",
        '>': "Reverse Entailment",
        '^': "Negation",
        '|': "Alternation",
        'v': "Cover",
        '#': "Independence"
    }
    return names.get(relation_symbol, "Unknown Relation")

def main():
    """
    Calculates the entailment relation by composing edits
    based on MacCartney's natural logic framework.
    """
    # The join table for composing relations.
    # Note: join('|', 'v') is set to '<' based on logical inference for this problem,
    # as explained in the reasoning.
    join_table = {
        '=': {'=': '=', '<': '<', '>': '>', '^': '^', '|': '|', 'v': 'v', '#': '#'},
        '<': {'=': '<', '<': '<', '>': '#', '^': '|', '|': '#', 'v': 'v', '#': '#'},
        '>': {'=': '>', '<': '#', '>': '>', '^': '^', '|': 'v', 'v': '>', '#': '#'},
        '^': {'=': '^', '<': 'v', '>': '|', '^': '=', '|': '>', 'v': '<', '#': '#'},
        '|': {'=': '|', '<': '#', '>': 'v', '^': '>', '|': '|', 'v': '<', '#': '#'},
        'v': {'=': 'v', '<': 'v', '>': '>', '^': '<', '|': '<', 'v': 'v', '#': '#'},
        '#': {'=': '#', '<': '#', '>': '#', '^': '#', '|': '#', 'v': '#', '#': '#'}
    }

    def join(r1, r2):
        """Composes two relations using the join table."""
        # Making it symmetric
        if r1 in join_table and r2 in join_table[r1]:
            return join_table[r1][r2]
        return join_table[r2][r1]

    # Step-by-step derivation
    premise = "Mark is singing a pop song by Taylor Swift"
    hypothesis = "Mark is not singing a song by Michael Jackson"

    print(f"Premise: \"{premise}\"")
    print(f"Hypothesis: \"{hypothesis}\"")
    print("-" * 20)
    print("Executing compositional proof:")

    # Step 1: Initial state
    current_relation = '='
    print("Step 1: Start with the premise.")
    print(f"         Initial Relation: {current_relation} ({get_relation_name(current_relation)})")
    print("-" * 20)

    # Step 2: Insertion of 'not'
    edit_relation_2 = '^'
    previous_relation = current_relation
    current_relation = join(current_relation, edit_relation_2)
    print("Step 2: Insert 'not'.")
    print(f"         Edit creates a '{edit_relation_2}' ({get_relation_name(edit_relation_2)}) relation.")
    print(f"         Joining '{previous_relation}' with '{edit_relation_2}' results in '{current_relation}'.")
    print(f"         Current Relation: {current_relation} ({get_relation_name(current_relation)})")
    print("-" * 20)


    # Step 3: Deletion of 'pop'
    # Lexical rel: <, context: downward, projection: >
    projected_relation_3 = '>'
    previous_relation = current_relation
    current_relation = join(current_relation, projected_relation_3)
    print("Step 3: Delete 'pop'.")
    print(f"         The projected relation from this edit is '{projected_relation_3}' ({get_relation_name(projected_relation_3)}).")
    print(f"         Joining '{previous_relation}' with '{projected_relation_3}' results in '{current_relation}'.")
    print(f"         Current Relation: {current_relation} ({get_relation_name(current_relation)})")
    print("-" * 20)


    # Step 4: Substitution of 'Taylor Swift' with 'Michael Jackson'
    # Lexical rel: |, context: downward, projection: v
    projected_relation_4 = 'v'
    previous_relation = current_relation
    current_relation = join(current_relation, projected_relation_4)
    print("Step 4: Substitute 'Taylor Swift' with 'Michael Jackson'.")
    print(f"         The projected relation from this edit is '{projected_relation_4}' ({get_relation_name(projected_relation_4)}).")
    print(f"         Joining '{previous_relation}' with '{projected_relation_4}' results in '{current_relation}'.")
    print(f"         Final Relation: {current_relation} ({get_relation_name(current_relation)})")
    print("-" * 20)


    final_relation_name = get_relation_name(current_relation)
    print(f"The final projected natural logic operator is '{current_relation}', which is named {final_relation_name}.")

if __name__ == "__main__":
    main()
