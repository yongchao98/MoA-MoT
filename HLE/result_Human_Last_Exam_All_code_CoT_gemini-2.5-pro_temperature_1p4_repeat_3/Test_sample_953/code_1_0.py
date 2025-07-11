def solve_mis_complexity():
    """
    Determines the complexity categories for the MIS algorithm on different graph classes.

    The analysis is based on established results in theoretical computer science:
    1.  For cycles and trees (which are planar graphs), the algorithm's runtime is O((log log n)^2).
        This complexity falls into category 6.
    2.  For general graphs with bounded degree, the runtime is Theta(log n).
        This complexity falls into category 9.
    """
    # d1: Complexity for a cycle of length n.
    # A cycle is a planar graph. The complexity is O((log log n)^2).
    # This falls into category 6.
    d1 = 6

    # d2: Complexity for a tree on n vertices with degree at most 100.
    # A tree is a planar graph. The complexity is O((log log n)^2).
    # This falls into category 6.
    d2 = 6

    # d3: Complexity for any graph on n vertices with degree at most 100.
    # For general graphs, the complexity is Theta(log n).
    # This falls into category 9.
    d3 = 9

    # The final result is the concatenation of these digits.
    final_answer_string = f"{d1}{d2}{d3}"
    print(final_answer_string)

solve_mis_complexity()