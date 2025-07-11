def solve():
    """
    Solves the problem by determining for each class of sets whether a maximal element always exists.
    The logic for each case is as follows:
    A: Yes. By Zorn's Lemma, the union of a chain of H-free graphs is H-free, providing an upper bound.
    B: Yes. Any non-empty finite set of real numbers has a maximum element.
    C: Depends. The set of natural numbers has no maximal element, but the set of negative integers does.
    D: Yes. The class of uncountable discrete subsets of R is empty, so the property is vacuously true.
    E: Yes. (a<=b if b is subseq of a). Constant sequences are maximal elements.
    F: No. (a<=b if a is subseq of b). For any sequence m, one can construct a sequence x of which m is a proper subsequence.
    """
    answer = "YYDYYN"
    print(answer)

solve()