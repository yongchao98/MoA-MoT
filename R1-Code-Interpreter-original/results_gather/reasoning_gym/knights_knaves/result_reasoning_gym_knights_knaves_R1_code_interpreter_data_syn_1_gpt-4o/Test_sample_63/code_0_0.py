def check_statements():
    # Scenario 1: Aurora is an angel
    aurora_is_angel = True
    owen_is_devil = not aurora_is_angel  # Based on Aurora's statement
    owen_statement = aurora_is_angel or not owen_is_devil  # Owen's statement
    if aurora_is_angel and not owen_statement:
        return "Contradiction in Scenario 1"

    # Scenario 2: Aurora is a devil
    aurora_is_angel = False
    owen_is_angel = not aurora_is_angel  # Based on Aurora's statement
    owen_statement = aurora_is_angel or owen_is_angel  # Owen's statement
    if not aurora_is_angel and owen_statement:
        return "Scenario 2 is consistent"

    return "No consistent scenario found"

print(check_statements())