def solve():
    """
    Solves the problem by determining for six classes of preordered sets
    whether they always have a maximal element.

    The analysis for each case is as follows:
    A) Depends on H. For H=K1, Yes. For H=K2, No. -> D
    B) Depends on the set. The empty set has no maximal element, but {1} does. -> D
    C) Depends on the set. Z has no maximal element, but {-n | n>=0} does. -> D
    D) The class of uncountable discrete subsets of R is empty, so the property is vacuously true. -> Y
    E) The relation is 'b is a subsequence of a'. Constant sequences are maximal. -> Y
    F) The relation is 'a is a subsequence of b'. No sequence is maximal. -> N
    """
    answer = "DDDYYN"
    print(answer)

solve()