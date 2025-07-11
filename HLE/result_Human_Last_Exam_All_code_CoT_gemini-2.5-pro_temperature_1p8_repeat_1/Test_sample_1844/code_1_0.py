def solve_paper_folding_problem():
    """
    Calculates the total number of edges on a piece of paper after being
    folded, cut, and unfolded according to the problem description.
    """

    # Step 1: Calculate the number of layers.
    # The paper is folded 4 times. Each fold doubles the layers.
    num_folds = 4
    num_layers = 2 ** num_folds
    print(f"The paper is folded 4 times, resulting in {num_layers} layers.")

    # Step 2: Analyze the cuts.
    # The four corners of the final folded square are cut.
    num_cuts = 4
    print(f"The {num_cuts} corners of the folded square are cut.")
    print("Because the folding is symmetric, all cuts occur in the interior of the original paper, creating holes.")
    
    # Step 3: Calculate the number of new edges from the cuts.
    # Each cut on the folded paper creates an edge on each of the layers it passes through.
    new_edges = num_cuts * num_layers
    print(f"Each of the {num_cuts} cuts creates an edge on all {num_layers} layers.")
    print(f"Number of new edges from the cuts = {num_cuts} * {num_layers} = {new_edges}.")

    # Step 4: Count the original edges.
    # Since the cuts are internal, the 4 edges of the original square are not cut.
    original_edges = 4
    print(f"The {original_edges} outer edges of the original square remain intact.")

    # Step 5: Calculate the total number of edges.
    # The total is the sum of the new edges from the holes and the original outer edges.
    total_edges = new_edges + original_edges
    print("\nThe total number of edges is the sum of the new edges and the original edges.")
    print(f"Final Equation: {new_edges} + {original_edges} = {total_edges}")

solve_paper_folding_problem()