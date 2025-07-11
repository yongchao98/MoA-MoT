import math

def count_h3_elements():
    """
    Calculates the number of elements of the reflection group of type H3
    that have a regular eigenvector with a corresponding eigenvalue of order 10.
    """
    # The group is H3, its Coxeter number h is 10.
    # We are looking for elements that are conjugate to c^k where c is a Coxeter element
    # and gcd(k, 10) = 1. The values for k are {1, 3, 7, 9}.

    # This requires us to find the size of the union of the conjugacy classes of
    # c^1, c^3, c^7, and c^9.

    # In H3, there are two distinct conjugacy classes of elements of order 10.
    # The Coxeter element c and its powers c^k (for k coprime to 10) are all
    # elements of order 10.

    # Let's call the two classes Class_A and Class_B.
    # A Coxeter element c belongs to one, say Class_A. Its inverse c^9 is also in Class_A.
    # The powers c^3 and its inverse c^7 belong to the other class, Class_B.

    # The size of a conjugacy class containing a Coxeter element c is |H3| / h.
    # |H3| = 120, h = 10.
    class_size = 120 / 10
    
    # Both conjugacy classes of order 10 elements in H3 have this size.
    size_class_A = int(class_size)
    size_class_B = int(class_size)
    
    # The elements we want are the union of these two disjoint classes.
    # So we sum their sizes.
    total_elements = size_class_A + size_class_B

    print("The total number of such elements is the sum of the sizes of two distinct conjugacy classes of order 10.")
    print(f"The first class, containing c^1 and c^9, has size {size_class_A}.")
    print(f"The second class, containing c^3 and c^7, has size {size_class_B}.")
    print("The final calculation is:")
    print(f"{size_class_A} + {size_class_B} = {total_elements}")

count_h3_elements()