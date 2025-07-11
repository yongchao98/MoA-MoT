from fractions import Fraction

def calculate_total_mass():
    """
    Calculates the total mass M(A_5, rho, 2) using Kedlaya's formula.
    """
    p = 2

    # Data for conjugacy classes of subgroups of A5.
    # Format: (name, N_H, H_order, homs_z2_h, element_data)
    # N_H: number of subgroups in the class
    # H_order: order of the subgroup
    # homs_z2_h: number of elements of 2-power order in H
    # element_data: list of (count, nfix) tuples for elements in H
    #   - count: number of elements of this type
    #   - nfix: number of fixed points for this permutation
    subgroup_classes = [
        ('Trivial', 1, 1, 1, [(1, 5)]),
        ('C2', 15, 2, 2, [(1, 5), (1, 1)]),
        ('C3', 10, 3, 1, [(1, 5), (2, 2)]),
        ('C5', 6, 5, 1, [(1, 5), (4, 0)]),
        ('V4', 5, 4, 4, [(1, 5), (3, 1)]),
        ('S3', 10, 6, 4, [(1, 5), (2, 2), (3, 1)]),
        ('D10', 6, 10, 6, [(1, 5), (4, 0), (5, 1)]),
        ('A4', 5, 12, 4, [(1, 5), (8, 2), (3, 1)]),
        ('A5', 1, 60, 16, [(1, 5), (20, 2), (15, 1), (24, 0)])
    ]

    total_mass = Fraction(0)
    terms_str = []

    for name, N_H, H_order, homs_z2_h, element_data in subgroup_classes:
        
        # Calculate the inner sum: sum_{g in H} p^(-codim(V^g))
        inner_sum = Fraction(0)
        for count, nfix in element_data:
            codim = 5 - nfix
            inner_sum += count * Fraction(1, p**codim)

        # Calculate the contribution of this subgroup class
        term = Fraction(N_H * homs_z2_h, H_order**2) * inner_sum
        total_mass += term
        terms_str.append(str(term))

    # Print the full equation
    print("The total mass is the sum of contributions from each subgroup class:")
    print(f"M(A5, rho, 2) = {' + '.join(terms_str)}")
    print("\nCalculating the sum:")
    print(f"M(A5, rho, 2) = {total_mass.numerator}/{total_mass.denominator}")

calculate_total_mass()
<<<52393/3600>>>