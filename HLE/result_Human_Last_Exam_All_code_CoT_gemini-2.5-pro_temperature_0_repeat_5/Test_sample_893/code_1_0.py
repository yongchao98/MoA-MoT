def solve():
    """
    Solves the problem by determining for each class of sets whether a maximal element always exists.
    A) H-free graphs: Yes, by Zorn's Lemma.
    B) Finite subset of R: Yes, has a maximum.
    C) Countable discrete subset of R: Depends (e.g., N vs Z-).
    D) Uncountable discrete subset of R: The class is empty, so the universal statement is vacuously true. Yes.
    E) a <= b if b is a subsequence of a: Yes, constant sequences are maximal.
    F) a <= b if a is a subsequence of b: No, for any sequence m, we can construct a strictly greater one.
    """
    answer = "YYDYYN"
    print(answer)

solve()