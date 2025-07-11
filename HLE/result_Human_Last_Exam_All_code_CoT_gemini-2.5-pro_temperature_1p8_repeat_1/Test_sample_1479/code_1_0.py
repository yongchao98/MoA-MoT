def explain_t2_vs_t3_structure():
    """
    Illustrates the structural difference between Type-2 and Type-3 fuzzy sets
    by showing the variables needed to define a membership grade.
    """
    # A single point in the primary universe of discourse
    primary_variable_x = 5.0
    
    print("This script demonstrates the dimensional difference between Type-2 and Type-3 fuzzy sets.")
    print(f"We will inspect the membership structure for a single input value: x = {primary_variable_x}\n")

    # --- Type-2 Fuzzy Set Structure ---
    print("--- Type-2 Fuzzy Membership Function ---")
    print("A Type-2 MF's domain has two variables: the primary variable (x) and the primary membership (u).")
    print("The function maps a point (x, u) to a secondary grade.")
    print("Equation form: secondary_grade = f(x, u)\n")

    # A discrete representation of the secondary MF for our specific x.
    # This is a "vertical slice" of the 3D Type-2 MF at x=5.0
    type2_membership_at_x = [
        {'u': 0.6, 'secondary_grade': 0.5},
        {'u': 0.7, 'secondary_grade': 1.0},
        {'u': 0.8, 'secondary_grade': 0.5}
    ]
    
    print(f"Example for x = {primary_variable_x}:")
    for point in type2_membership_at_x:
        u_val = point['u']
        grade = point['secondary_grade']
        print(f"f(x={primary_variable_x}, u={u_val}) = {grade}")
        
    print("-" * 30 + "\n")

    # --- Type-3 Fuzzy Set Structure ---
    print("--- Type-3 Fuzzy Membership Function ---")
    print("A Type-3 MF adds another layer. Its domain has THREE variables: x, u, and the secondary membership (v).")
    print("The function maps a point (x, u, v) to a tertiary grade.")
    print("Equation form: tertiary_grade = g(x, u, v)\n")

    # For a Type-3 set, the secondary grade from the Type-2 example is now a fuzzy set.
    # We will expand on the u=0.7 case from above. Its grade was 1.0. Now, that grade is fuzzy.
    u_of_interest = 0.7
    type3_membership_at_x_u = [
        {'v': 0.9, 'tertiary_grade': 0.8},
        {'v': 1.0, 'tertiary_grade': 1.0}
    ]

    print(f"Example for x = {primary_variable_x} and u = {u_of_interest}:")
    print(f"The secondary grade is no longer a single number, but a fuzzy set defined over variable 'v'.")
    for point in type3_membership_at_x_u:
        v_val = point['v']
        grade = point['tertiary_grade']
        print(f"g(x={primary_variable_x}, u={u_of_interest}, v={v_val}) = {grade}")
        
    print("\nConclusion: The fundamental structural change from Type-2 to Type-3 is the expansion from a two-variable domain (x, u) to a three-variable domain (x, u, v).")

# Execute the explanation function
explain_t2_vs_t3_structure()