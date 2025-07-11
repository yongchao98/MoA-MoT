def calculate_junction_chern_number():
    """
    Calculates the Chern number of a junction formed by two Chern insulators.
    """
    # Chern number of the first insulator
    c1 = 1
    
    # Chern number of the second insulator
    c2 = 1
    
    # The Chern number of the combined system is the sum of the individual Chern numbers.
    c_total = c1 + c2
    
    print("The Chern number of a junction of two topological insulators is the sum of their individual Chern numbers.")
    print(f"Chern number of the first insulator (C1): {c1}")
    print(f"Chern number of the second insulator (C2): {c2}")
    print("The total Chern number (C_total) of the junction is calculated as:")
    print(f"C_total = C1 + C2")
    # Outputting each number in the final equation as requested
    print(f"C_total = {c1} + {c2}")
    print(f"Therefore, the final Chern number is: {c_total}")

calculate_junction_chern_number()