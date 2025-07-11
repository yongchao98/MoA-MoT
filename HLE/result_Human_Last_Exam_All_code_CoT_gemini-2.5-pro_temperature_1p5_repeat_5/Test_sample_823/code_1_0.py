import sys

def explain_answer():
    """
    This function provides a step-by-step explanation for the graph theory problem.
    """
    print("Analyzing the properties of a class of graphs C with bounded degree and unbounded treewidth.")
    print("The task is to determine which statement must be true for such a class C.")
    print("Let's evaluate each option:\n")

    print("--- Option A: Contains induced cycle of length k ---")
    print("This is FALSE. The class of square grids has bounded degree (4) and unbounded treewidth, but only contains induced cycles of length 4.\n")

    print("--- Option B: Contains k-by-k grid as a subgraph ---")
    print("This is FALSE. A class of expander graphs has bounded degree and unbounded treewidth, but their structure prevents them from containing large grid subgraphs.\n")

    print("--- Option C: C is a family of expander graphs ---")
    print("This is FALSE. The class of square grids serves as a counterexample, as it meets the conditions but grids are not expanders.\n")

    print("--- Option D: Contains induced matching of size k ---")
    print("This is not guaranteed. Whether large treewidth implies a large induced matching is an open research problem. Therefore, it cannot be something that *must* be true.\n")

    print("--- Option E: Contains clique-minors of unbounded size ---")
    print("This is TRUE. This follows from fundamental theorems in structural graph theory.\n")
    print("The logical implication is as follows:")

    # There are no numbers in an equation, but we can represent the logical steps.
    print("1. Unbounded Treewidth in C => For any r, exists G in C with an r x r grid minor. (by the Excluded Grid Theorem)")
    print("2. An r x r grid graph => Contains a K_k clique minor where k is proportional to r (e.g., k = floor(r/2)).")
    print("3. By transitivity of the minor relation: a graph in C with a large grid minor must also have a large clique minor.")
    
    print("\nSince the treewidth in C is unbounded, the size of the grid minors is unbounded, and therefore the size of the clique minors is also unbounded.")
    print("Conclusion: Statement E must be true.")

if __name__ == "__main__":
    explain_answer()