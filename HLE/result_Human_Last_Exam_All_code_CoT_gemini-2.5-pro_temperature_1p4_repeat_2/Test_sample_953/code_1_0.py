import math

def solve():
    """
    Analyzes the runtime of a parallel MIS algorithm on three graph classes.

    The algorithm is a variant of Luby's algorithm. For graphs with a maximum degree
    bounded by a constant (e.g., 100), it's a standard result in parallel
    algorithms that this type of algorithm terminates in O(log n) steps with
    high probability. This provides an upper bound for all three cases.

    To find the tight complexity, we also need a lower bound. We can find a "worst-case"
    graph in each class that forces the algorithm to take Omega(log n) steps.

    1. f_1(n) for a cycle of length n:
       - A cycle is a 2-regular graph (max degree is 2).
       - Let's analyze the probability that a vertex `v` survives one step. A vertex
         `v` survives if it's not chosen and none of its neighbors are chosen. For
         a cycle, a vertex `v` has two neighbors, `u` and `w`.
       - The probability of `v` surviving can be calculated to be a constant `c < 1`
         (specifically, 2/15).
       - This means that in each step, the number of remaining vertices is expected
         to decrease by a constant factor. Such a process takes O(log n) steps
         to complete.
       - A more rigorous analysis shows this bound is tight. Therefore,
         f_1(n) = Theta(log n).

    2. f_2(n) for any tree on n vertices of degree at most 100:
       - This class of graphs includes a simple path, which is a tree with max degree 2.
       - The analysis for a long path is identical to that for a cycle for any vertex
         not near an endpoint.
       - Since a path requires Theta(log n) steps, the function f_2(n), which must
         be an upper bound for *all* trees, must be at least Omega(log n).
       - Combined with the general O(log n) upper bound for bounded-degree graphs,
         we have f_2(n) = Theta(log n).

    3. f_3(n) for any graph on n vertices of degree at most 100:
       - This class is the most general. It includes cycles and paths as special cases.
       - Since a cycle C_n is a graph with max degree 2 (and thus <= 100), the
         argument for f_1(n) shows that f_3(n) must be at least Omega(log n).
       - The general O(log n) upper bound holds for this class.
       - Therefore, f_3(n) = Theta(log n).

    Mapping to categories:
    The function f(n) = Theta(log n) does not fit in categories 1 through 8.
    For example, category 8 requires f(n) = o(log n), which is not true for
    f(n) = Theta(log n).
    Category 9 is f(n) = Omega(log n). Since Theta(log n) implies Omega(log n),
    this is the correct category for all three cases.

    - d_1 (for cycles) = 9
    - d_2 (for trees) = 9
    - d_3 (for general graphs) = 9
    """

    # For f_1(n) on a cycle, the runtime is Theta(log n). This is Omega(log n).
    d_1 = 9
    # For f_2(n) on a bounded-degree tree, the runtime is Theta(log n). This is Omega(log n).
    d_2 = 9
    # For f_3(n) on a bounded-degree graph, the runtime is Theta(log n). This is Omega(log n).
    d_3 = 9

    # The final answer is the concatenation of the three digits.
    final_answer_string = f"{d_1}{d_2}{d_3}"
    print(final_answer_string)

solve()
<<<999>>>