def solve():
    """
    This function solves the mathematical problem posed by the user.

    The problem asks for the number of distinct homeomorphism classes for a compact
    topological space X with the following properties:
    1. X contains a dense copy of the long ray R.
    2. Every bounded continuous function f: R -> R extends to a unique continuous
       function on X.

    These two properties are the defining characteristics of the Stone-Čech
    compactification of R, denoted as βR.

    The Stone-Čech compactification is unique up to homeomorphism. This means that
    any two spaces that satisfy these properties must be homeomorphic to each other.

    Therefore, all such spaces X belong to the same homeomorphism class.
    The number of distinct homeomorphism classes is 1.
    """
    # The "equation" is simply that the result is 1.
    number_of_classes = 1
    
    # Output the single number from the "equation"
    print(number_of_classes)

solve()