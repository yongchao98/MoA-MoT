def solve():
    """
    This function determines the answer for each case and prints the final combined string.
    A: N - For any H-free graph, adding a disjoint vertex creates a larger H-free graph.
    B: D - A non-empty finite set of reals has a max, but the empty set doesn't.
    C: D - {1,2,3} has a max, but N = {1,2,3,...} doesn't.
    D: Y - The class of uncountable discrete subsets of R is empty, so the condition is vacuously true.
    E: Y - Constant sequences (c, c, c, ...) are maximal elements.
    F: N - For any sequence m, one can construct a strictly greater sequence x.
    """
    answer = "NDDYYN"
    print(answer)

solve()