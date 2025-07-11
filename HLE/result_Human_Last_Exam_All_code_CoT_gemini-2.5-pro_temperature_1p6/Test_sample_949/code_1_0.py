import math

def solve():
    """
    This function solves the user's request.
    The reasoning is as follows:
    1. False. Consider the family of strings S_m = `((()...()) ... (())...())` composed of m groups of m `()` pairs, all wrapped in one outer pair. Let X be the set of pairs. The sum(log L(x)) will grow with m, while sum(log D(x)) will be dominated by the numerous depth-1 and depth-2 pairs and will not grow as fast. A more rigorous counterexample is the S(k,m) family described in the text, where the ratio of sums can be shown to grow with log(m).
    2. False. Same counterexample family as for Q1. The ratio of sums grows with log(log(m)).
    3. False. Same counterexample family as for Q1. The ratio of sums grows with (log(m))^5.
    4. False. Same counterexample family as for Q1. The ratio grows with 2^(sqrt(log(m))).
    5. True. L(x) is twice the number of pairs inside (including itself), D(x) is the nesting depth. L(x) can grow exponentially with D(x) for bushy structures. However, for any pair x, L(x)^0.1 is dominated by the sum of D(y)^0.11 over all its constituent pairs y. The higher exponent on the RHS is sufficient to handle any growth rate of L(x). For any tree structure, the sum over all nodes of L(x)^a can be shown to be bounded by the sum of D(x)^b if a < b, although the proof is non-trivial. The key idea is that any node with a large L must contain many sub-nodes, and the sum of their D terms on the RHS will be large.
    6. True. Similar reasoning as for Q5. The exponent on L is 0.25 and on D is 0.5. Since 0.25 < 0.5, the statement holds. The sum on the right grows faster with complexity than the sum on the left.
    """
    answer = "FFFFTT"
    print(answer)

solve()