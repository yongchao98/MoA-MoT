def solve_enclitics_order():
    """
    This function determines and prints the correct order of a given set of
    Old Russian enclitics based on established linguistic rules.
    """
    # The given enclitics: бы, же, еси, бо, мя

    # The correct order is determined by linguistic rules for enclitic clusters in Old Russian:
    # 1. Modal Particle (бы)
    # 2. Discourse Particles (же, бо - in that order)
    # 3. Pronominal Accusative (мя)
    # 4. Verbal (еси)
    correct_order = ['бы', 'же', 'бо', 'мя', 'еси']

    # Prepare the output string
    result_string = ", ".join(correct_order)
    
    print("The correct order of the enclitics is:")
    print(result_string)

solve_enclitics_order()
<<<бы, же, бо, мя, еси>>>