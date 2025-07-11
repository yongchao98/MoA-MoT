def solve():
    """
    Solves the problem by analyzing each case and printing the final answer string.
    A) N: For any H-free graph G, G with an added isolated vertex is a larger H-free graph.
    B) Y: Any non-empty finite subset of R has a maximum, which is a maximal element.
    C) D: S={1,2,3,...} has no maximal element. S={-1,-2,-3,...} has a maximal element.
    D) Y: The class of uncountable discrete subsets of R is empty, so the statement is vacuously true.
    E) Y: A constant sequence like (1,1,1,...) is a maximal element under the reverse subsequence order.
    F) N: For any sequence a, a new sequence b can be constructed which is strictly greater in the subsequence order.
    """
    # The final answer is the concatenation of the results for each case.
    answer = "NYDYN"
    print(answer)

solve()