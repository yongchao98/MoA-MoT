import fractions

def demonstrate_topology_classes():
    """
    This function demonstrates the concepts by printing examples of subsets of Q.
    """
    # Part 1: Illustrate two sets that are homeomorphic to a subset of the other.
    # We choose A = Q intersect [0, 1] and B = Q intersect [0, 2].
    # A is a subset of B. The embedding is f(x) = x.
    # B can be mapped to a subset of A via g(x) = x / 2.
    # So, by the topological Schroeder-Bernstein theorem, they are homeomorphic.
    
    print("Part 1: Example of Two Equivalent Subsets")
    print("Let A = Q intersect [0, 1] and B = Q intersect [0, 2].")
    print("1. An embedding f: A -> B is given by f(x) = x.")
    # Example:
    a_point = fractions.Fraction(1, 2)
    print(f"   f({a_point}) = {a_point}, which is in B.")
    print("2. An embedding g: B -> A is given by g(x) = x / 2.")
    b_point = fractions.Fraction(3, 2)
    print(f"   g({b_point}) = {b_point / 2}, which is in A.")
    print("Conclusion: A and B belong to the same equivalence class (are homeomorphic).\n")

    # Part 2 & 3: Illustrating different equivalence classes (homeomorphism types)
    # The number of such classes is 2^{\aleph_0}, the cardinality of the continuum.
    # We illustrate some of the simpler classes by generating their members.
    
    print("Part 2 & 3: Examples of Different Equivalence Classes\n")

    def print_set(name, description, generator_func, limit=15):
        """Helper function to print examples of sets."""
        print(f"Class Type: {name}")
        print(f"Description: {description}")
        # Generate and print each number in the set explicitly
        points = generator_func(limit)
        points_str = [str(p) for p in points]
        print(f"Example points: {{ {', '.join(points_str)}, ... }}")
        print("-" * 40)

    # Type 1: Finite sets. One class for each size n.
    def gen_finite_set(n):
        return [fractions.Fraction(i) for i in range(1, n + 1)]

    print_set("Finite (n=5)", "A finite set with 5 points.",
              lambda limit: gen_finite_set(5))

    # Type 2: Countably infinite discrete set.
    def gen_discrete_infinite(n):
        return [fractions.Fraction(i) for i in range(n)]

    print_set("Discrete Infinite", "Homeomorphic to the integers N.",
              gen_discrete_infinite)

    # Type 3: A convergent sequence. (Cantor-Bendixson rank 2)
    def gen_convergent_sequence(n):
        return [fractions.Fraction(1, k) for k in range(1, n + 1)] + [fractions.Fraction(0)]

    print_set("Convergent Sequence", "A sequence {1/k} and its limit point 0.",
              gen_convergent_sequence)

    # Type 4: Disjoint union of k convergent sequences.
    def gen_multi_convergent(k, points_per_seq):
        result = []
        for i in range(k):
            limit_point = fractions.Fraction(i)
            # Create a sequence converging to i, e.g., {i + 1/j}
            seq = [limit_point + fractions.Fraction(1, j) for j in range(1, points_per_seq + 1)]
            result.extend(seq)
            result.append(limit_point)
        return result

    print_set("Two Convergent Sequences", "Two disjoint copies of a convergent sequence.",
              lambda limit: gen_multi_convergent(2, (limit-2) // 2))

    # Type 5: A space of higher rank (Cantor-Bendixson rank 3)
    def gen_rank_3_space(n_limit_pts):
        """Generates a space S where S' is a convergent sequence."""
        result = []
        # The set of limit points S' is {1/k} U {0}
        level2_limit_pts = [fractions.Fraction(1, k) for k in range(1, n_limit_pts + 1)]
        level3_limit_pt = fractions.Fraction(0)
        
        # Add S' to the set S
        result.extend(level2_limit_pts)
        result.append(level3_limit_pt)

        # For each point in S' (except 0), add a sequence converging to it
        for p in level2_limit_pts:
            # Add points p + 1/m where m is large enough to avoid collisions
            for j in range(1, 4):
                 result.append(p + fractions.Fraction(1, p.denominator * 1000 + j))
        return result

    print_set("Convergent Sequence of Sequences", "Points converging to points of another convergent sequence.",
              lambda limit: gen_rank_3_space(limit // 4))

    print("And so on. The variety of such structures is vast.")
    print("The total number of non-homeomorphic classes is the cardinality of the continuum.")

demonstrate_topology_classes()