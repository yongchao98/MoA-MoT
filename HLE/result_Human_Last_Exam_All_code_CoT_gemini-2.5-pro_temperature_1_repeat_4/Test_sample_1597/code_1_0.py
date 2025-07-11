def count_allowed_pairs():
    """
    Calculates the number of "allowed" ordered pairs (a,b) with 1 <= a, b <= 1000.

    An ordered pair (a,b) is allowed if for every primitive function f,
    f(ab) >= f(a)f(b). This condition is equivalent to a|b or b|a.
    """
    limit = 1000

    # Calculate |S1|, the number of pairs (a,b) where a|b.
    # For each 'a', the number of 'b's that are multiples of 'a' up to the limit
    # is floor(limit/a).
    s1_count = 0
    for a in range(1, limit + 1):
        s1_count += limit // a

    # |S2| = |S1| by symmetry (number of pairs where b|a).
    s2_count = s1_count

    # |S1 intersect S2| is the number of pairs where a|b and b|a, which means a=b.
    # These are (1,1), (2,2), ..., (1000,1000).
    intersection_count = limit

    # Using the Principle of Inclusion-Exclusion: |S1 U S2| = |S1| + |S2| - |S1 intersect S2|
    total_allowed_pairs = s1_count + s2_count - intersection_count

    # Output the components of the final calculation.
    print(f"2 * {s1_count} - {intersection_count} = {total_allowed_pairs}")

count_allowed_pairs()
