def check_knight_knave():
    # Scenario 1: Ava is a knight
    ava_knight = True
    scarlett_knave = not ava_knight
    scarlett_statement = ava_knight  # If Ava is a knight, Scarlett is a knight
    if ava_knight and scarlett_knave and not scarlett_statement:
        return "Ava is a knight, and Scarlett is a knave"

    # Scenario 2: Ava is a knave
    ava_knight = False
    scarlett_knave = not ava_knight
    scarlett_statement = ava_knight  # If Ava is a knight, Scarlett is a knight
    if not ava_knight and not scarlett_knave and not scarlett_statement:
        return "Ava is a knave, and Scarlett is a knight"

    return "No consistent scenario found"

print(check_knight_knave())