import sys

def solve():
    """
    This function solves for the minimum number of bipartite graphs required to cover the complete graph K_n.

    The problem asks for the "minimum number of vertices in a family of bipartite graphs"
    covering K_n. This is likely a slightly ambiguous phrasing for a classic problem in
    combinatorics. The standard problem asks for the minimum number of bipartite *graphs*
    in such a family.

    A literal interpretation ("minimum number of vertices in the union of vertex sets") would
    trivially result in n, since all vertices must be present to cover all edges incident to them.

    The more profound and likely intended question refers to the number of graphs. According to a
    theorem by Graham, Pollak, and Tverberg, the minimum number of bipartite graphs needed to
    cover all edges of a complete graph K_n is exactly n-1.

    The question asks to provide the final equation. For n=35, the equation is 35 - 1 = 34.
    """
    # The number of vertices in the complete graph K_n.
    n = 35

    # According to the theorem, the minimum number of bipartite graphs required is n - 1.
    result = n - 1

    # Print the final equation as requested.
    print(f"{n} - 1 = {result}")

solve()