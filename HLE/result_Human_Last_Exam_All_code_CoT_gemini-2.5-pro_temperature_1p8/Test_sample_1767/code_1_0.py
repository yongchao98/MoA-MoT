def calculate_junction_chern_number():
    """
    Calculates the Chern number of a junction of two Chern insulators.

    The Chern number is an additive topological invariant. When two topological
    insulators are joined, the Chern number of the composite system is the sum
    of their individual Chern numbers.
    """
    # Chern number of the first insulator
    c1 = 1
    # Chern number of the second insulator
    c2 = 1

    # The total Chern number of the composite system is the sum of the individuals.
    c_total = c1 + c2

    print(f"The system consists of two Chern insulators joined together.")
    print(f"The first insulator has a Chern number C1 = {c1}.")
    print(f"The second insulator has a Chern number C2 = {c2}.")
    print("\nThe Chern number of the composite system is the sum of the individual Chern numbers:")
    print(f"C_total = C1 + C2")
    print(f"C_total = {c1} + {c2} = {c_total}")

calculate_junction_chern_number()
<<<2>>>