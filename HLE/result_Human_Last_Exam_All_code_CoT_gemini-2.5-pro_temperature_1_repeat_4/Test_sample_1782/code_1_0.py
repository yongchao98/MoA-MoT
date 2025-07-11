def diagonalization_demonstration():
    """
    This function illustrates the diagonalization argument used to prove
    the existence of the tree described in the problem.
    """
    print("This program demonstrates a key proof technique: Diagonalization.\n")

    # In the actual proof, this would be a list of all possible 'common refinements',
    # which has a cardinality of 2^omega_1. We can enumerate omega_1 of them.
    # Here, we represent them as a finite list of "forbidden paths".
    # A path is a sequence of choices (0 or 1 for a binary split).
    forbidden_paths = [
        [0, 1, 1, 0, 1, 0],
        [1, 0, 1, 1, 0, 0],
        [0, 0, 0, 1, 0, 1],
        [1, 1, 0, 0, 1, 1],
        [0, 1, 0, 1, 0, 1],
        [1, 0, 1, 0, 1, 0]
    ]
    
    num_levels = len(forbidden_paths)
    print(f"We have {num_levels} 'forbidden paths' we want to avoid:")
    for i, path in enumerate(forbidden_paths):
        print(f"  Path {i}: {path}")
    print("-" * 30)

    # We will construct a new path that is different from all forbidden paths.
    constructed_path = []

    print("Constructing a new path by diagonalizing:")
    for i in range(num_levels):
        # At level `i`, we look at the `i`-th choice of the `i`-th forbidden path.
        diagonal_element = forbidden_paths[i][i]
        
        # We choose the opposite of that choice for our new path.
        # This is the core of the diagonalization.
        # This choice ensures our path differs from forbidden_path[i] at position i.
        our_choice = 1 - diagonal_element
        
        # The 'equation' for our choice at this step
        print(f"Step {i}: The forbidden choice on the diagonal is forbidden_paths[{i}][{i}] = {diagonal_element}.")
        print(f"         Our choice is 1 - {diagonal_element} = {our_choice}")
        
        constructed_path.append(our_choice)

    print("-" * 30)
    print("The final constructed path is:")
    print(constructed_path)
    print("\nThis path is guaranteed to be different from every path in the forbidden list.")

    # Verification step
    for i in range(num_levels):
        print(f"Comparing with forbidden path {i}: It differs at index {i} ({constructed_path[i]} vs {forbidden_paths[i][i]})")

diagonalization_demonstration()
