def check_knight_knave():
    # Scenario 1: Owen is a knight, Charlotte is a knave
    owen_knight = True
    charlotte_knave = True
    owen_statement = charlotte_knave  # Owen says Charlotte is a knave
    charlotte_statement = not (not charlotte_knave or owen_knight)  # Charlotte's conditional statement

    scenario_1_valid = owen_knight == owen_statement and not charlotte_knave == charlotte_statement

    # Scenario 2: Owen is a knave, Charlotte is a knight
    owen_knight = False
    charlotte_knave = False
    owen_statement = charlotte_knave  # Owen says Charlotte is a knave
    charlotte_statement = not (not charlotte_knave or owen_knight)  # Charlotte's conditional statement

    scenario_2_valid = not owen_knight == owen_statement and not charlotte_knave == charlotte_statement

    if scenario_1_valid:
        return "Owen is a knight, and Charlotte is a knave"
    elif scenario_2_valid:
        return "Owen is a knave, and Charlotte is a knight"
    else:
        return "No valid scenario"

print(check_knight_knave())