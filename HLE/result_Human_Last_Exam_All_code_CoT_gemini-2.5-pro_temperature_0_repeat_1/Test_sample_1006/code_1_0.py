def solve():
    """
    This function calculates the number of distinct homeomorphism classes for the space X.

    The problem describes a compact space X which is a compactification of the long ray R.
    The key property is that every bounded continuous function from R to the real numbers
    extends uniquely to a continuous function on X. This is the defining universal property
    of the Stone-Čech compactification, denoted as βR.

    The Stone-Čech compactification of a Tychonoff space is unique up to homeomorphism.
    Therefore, any space X satisfying the given properties must be homeomorphic to βR.
    This implies that there is only one such homeomorphism class.
    """
    
    # The number of distinct homeomorphism classes is 1.
    num_classes = 1
    
    # The problem asks for the number of classes.
    # There is no equation, so we just print the final number.
    print(num_classes)

solve()