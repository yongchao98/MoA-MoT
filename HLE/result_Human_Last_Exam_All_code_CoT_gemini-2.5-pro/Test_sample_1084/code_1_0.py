import collections

def solve():
    """
    Determines the truthfulness of six statements about graph properties.

    The statements relate the number of 4-cycles (A) and specific 6-cycles (B)
    in a d-regular graph.

    Based on known results in extremal graph theory, the number of 4-cycles (A)
    and 6-cycles (B) are related. A general form of this relationship is:
    A <= c1 * (B/d) + c2 * n * d^2
    where c1 and c2 are constants. This means A is bounded by a term proportional
    to B/d and an error term. The error term nd^2 is a known upper bound for A
    in some families of graphs where B can be small or zero.

    Let's evaluate each statement based on the relationship A = O(B/d + nd^2).

    1) A = O(B/d^0.5 + nd^2):
       Since d >= 1, B/d <= B/d^0.5. The error term nd^2 is the same.
       So, O(B/d + nd^2) is a subset of O(B/d^0.5 + nd^2).
       This statement is TRUE.

    2) A = O(B/d^0.99 + nd^2):
       Since d >= 1, B/d <= B/d^0.99.
       This statement is TRUE.

    3) A = O(B/d + nd^2):
       This is the base relationship we are using.
       This statement is TRUE.

    4) A = O(B/d^1.5 + nd^2):
       This statement would imply that A grows slower than B/d.
       However, in random regular graphs, the number of 6-cycles is roughly d times
       the number of 4-cycles, meaning A is proportional to B/d.
       A counterexample would be a graph family where A is proportional to B/d,
       making the statement A = O(A/d^0.5 + nd^2), which is false for large d.
       This statement is FALSE.

    5) A = O(B/d^1.3 + n*d^2.6):
       Similar to statement 4, the exponent on d (1.3) is greater than 1.
       This statement is FALSE.

    6) A = O(B/d + n*d^1.5):
       The error term nd^1.5 is smaller than the known tight bound of nd^2.
       There exist graphs with A proportional to nd^2 but with very few 6-cycles (small B).
       For these graphs, the statement would become nd^2 = O(nd^1.5), which is false.
       This statement is FALSE.

    So, the results are Y, Y, Y, N, N, N.
    """
    results = ['Y', 'Y', 'Y', 'N', 'N', 'N']
    print("".join(results))

solve()