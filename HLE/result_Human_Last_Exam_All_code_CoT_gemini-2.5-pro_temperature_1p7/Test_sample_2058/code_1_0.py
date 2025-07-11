def solve():
    """
    Calculates the total mass M(A_5, rho, 2).
    """

    # Step 3: Calculate n(A_5, 2), the number of homomorphisms from Gamma_2 to A_5.
    # Data for subgroup classes of A_5:
    # H: (Number of subgroups N_H, Number of H-extensions of Q_2 a_H)
    subgroup_data = {
        'C1': (1, 1),
        'C2': (15, 7),
        'C3': (10, 4),
        'V4': (5, 19),
        'C5': (6, 1),
        'S3': (10, 9),
        'D5': (6, 1),
        'A4': (5, 1),
        'A5': (1, 1),
    }

    n_A5_2 = 0
    calculation_str = []
    for name, (N_H, a_H) in subgroup_data.items():
        term = N_H * a_H
        n_A5_2 += term
        calculation_str.append(f"({N_H} * {a_H})")

    # Order of A_5
    abs_A5 = 60

    # Final result is n(A_5, 2) / |A_5|
    # M(A_5, rho, 2) = n_A5_2 / abs_A5
    
    numerator = n_A5_2
    denominator = abs_A5

    # We need to express the result as a fraction in lowest terms.
    # Python's math.gcd can find the greatest common divisor.
    import math
    common_divisor = math.gcd(numerator, denominator)
    
    lowest_num = numerator // common_divisor
    lowest_den = denominator // common_divisor
    
    print(f"The number of subgroups and extensions for each class [H] are:")
    print(f"H=C1: N=1, a=1")
    print(f"H=C2: N=15, a=7")
    print(f"H=C3: N=10, a=4")
    print(f"H=V4: N=5, a=19")
    print(f"H=C5: N=6, a=1")
    print(f"H=S3: N=10, a=9")
    print(f"H=D5: N=6, a=1")
    print(f"H=A4: N=5, a=1")
    print(f"H=A5: N=1, a=1")
    print(f"\nThe total number of homomorphisms n(A_5, 2) is the sum of (N_H * a_H) over all classes [H]:")
    print(f"n(A_5, 2) = {' + '.join(calculation_str)}")
    print(f"n(A_5, 2) = {1*1} + {15*7} + {10*4} + {5*19} + {6*1} + {10*9} + {6*1} + {5*1} + {1*1}")
    print(f"n(A_5, 2) = {numerator}")
    print(f"\nThe order of the group A_5 is |A_5| = {abs_A5}.")
    print(f"The total mass M(A_5, rho, 2) is n(A_5, 2) / |A_5|.")
    print(f"M(A_5, rho, 2) = {numerator}/{denominator}")
    print(f"The fraction in lowest terms is {lowest_num}/{lowest_den}.")

solve()