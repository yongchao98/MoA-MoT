import textwrap

def solve_topology_problem():
    """
    Solves the problem by identifying the topological space X and using the uniqueness theorem.
    """

    explanation = """
    Let's analyze the properties of the compact topological space X.

    1.  The Problem Statement:
        - X is a compact topological space.
        - (1) X contains a dense copy of the long ray R = [0, \omega_1). This means X is a compactification of R.
        - (2) Every bounded continuous function f: R -> R extends to a unique continuous function on X.

    2.  Identifying the Space:
        The two properties given are the defining properties of the Stone-Čech compactification of the long ray R, which is denoted as \beta R.
        - Property (1) establishes that X is a compactification of R.
        - Property (2) is the universal mapping property that uniquely characterizes the Stone-Čech compactification among all other compactifications of R.

    3.  Uniqueness of the Stone-Čech Compactification:
        A fundamental theorem in topology states that the Stone-Čech compactification of a Tychonoff space (like the long ray R) is unique up to homeomorphism. This means that any two spaces that satisfy the given conditions must be topologically identical (homeomorphic).

    4.  Conclusion:
        Since any space X that satisfies the given properties must be homeomorphic to the Stone-Čech compactification of R (\beta R), all such spaces belong to the same homeomorphism class. Therefore, there is only one distinct homeomorphism class.

    5.  Calculation:
        The number of distinct homeomorphism classes is determined by this uniqueness.
    """

    print(textwrap.dedent(explanation).strip())

    num_classes = 1

    print("\nFinal Calculation:")
    print(f"Number of classes = {num_classes}")


solve_topology_problem()

print("\n<<<1>>>")