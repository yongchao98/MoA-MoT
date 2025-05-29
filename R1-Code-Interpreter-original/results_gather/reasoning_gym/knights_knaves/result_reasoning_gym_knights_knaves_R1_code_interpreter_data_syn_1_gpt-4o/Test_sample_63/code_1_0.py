def check_statements():
    # Case 1: Assume Aurora is an angel
    aurora_is_angel = True
    owen_is_devil = not aurora_is_angel  # From Aurora's statement
    owen_statement = aurora_is_angel or not owen_is_devil  # Owen's statement
    if owen_is_devil and not owen_statement:
        return "Aurora is an angel, and Owen is a devil"

    # Case 2: Assume Aurora is a devil
    aurora_is_angel = False
    owen_is_angel = not aurora_is_angel  # From Aurora's statement
    owen_statement = aurora_is_angel or owen_is_angel  # Owen's statement
    if owen_is_angel and owen_statement:
        return "Aurora is a devil, and Owen is an angel"

    return "No consistent solution found"

print(check_statements())