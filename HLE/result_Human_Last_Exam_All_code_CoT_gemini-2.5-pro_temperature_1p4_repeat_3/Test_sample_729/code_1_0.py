def solve_power_subgroups():
    """
    Calculates and explains the number of power subgroups in the
    generalized quaternion group of size 128.
    """
    group_size = 128
    n = group_size // 4

    print(f"Let G be the generalized quaternion group of size {group_size}, denoted Q_{group_size}.")
    print(f"The presentation of this group is <x, y | x^{2*n} = 1, x^n = y^2, y^(-1)xy = x^(-1)>.")
    print(f"For Q_{group_size}, we have n = {n}, so the presentation is <x, y | x^64 = 1, x^32 = y^2, y^(-1)xy = x^(-1)>.\n")

    print("A power subgroup G^k is the set of all k-th powers of elements in G.")
    print("To find all unique power subgroups, we can analyze G^d for each divisor d of the group's exponent.\n")

    print("Step 1: Find the exponent of G = Q_128.")
    print("The elements are of the form x^i (orders divide 64) and x^i*y (order 4).")
    exponent = 64
    print(f"The exponent is the least common multiple of all element orders, which is lcm(64, 4) = {exponent}.\n")

    print(f"Step 2: Find all divisors of the exponent, {exponent}.")
    divisors = sorted([d for d in range(1, exponent + 1) if exponent % d == 0])
    print(f"The divisors of {exponent} are: {divisors}.\n")

    print("Step 3: Characterize the power subgroup G^d for each divisor d.")
    subgroups = []
    for d in divisors:
        if d == 1:
            name = f"G = Q_{group_size}"
            order = group_size
        else:
            name = f"<x^{d}>"
            order = (2 * n) // d
        subgroups.append({'d': d, 'name': name, 'order': order})

    print("The distinct power subgroups are:")
    for sub in subgroups:
        print(f"- G^{sub['d']} is the subgroup {sub['name']}, with order {sub['order']}.")
    
    count = len(subgroups)
    print(f"\nAll {count} of these subgroups have distinct orders, so they are unique.")
    
    equation_parts = ["1"] * count
    final_equation = " + ".join(equation_parts)
    
    print("\nFinal calculation:")
    print(f"Counting each unique subgroup: {final_equation} = {count}")


solve_power_subgroups()
