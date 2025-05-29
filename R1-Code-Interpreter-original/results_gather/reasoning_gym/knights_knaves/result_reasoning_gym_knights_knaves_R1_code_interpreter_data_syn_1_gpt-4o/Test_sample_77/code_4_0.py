def check_statements():
    # Assume Victoria is a laggard (liar)
    victoria_pioneer = False
    amelia_pioneer = False  # Assume Amelia is a laggard (liar)

    # Victoria's statement: If Amelia is a pioneer, then Amelia is a laggard
    victoria_statement = not amelia_pioneer or not amelia_pioneer

    # Amelia's statement: Victoria is a pioneer if and only if Amelia is a pioneer
    amelia_statement = (victoria_pioneer == amelia_pioneer)

    # Check if both statements are consistent with the assumptions
    if not victoria_statement and not amelia_statement:
        return "Victoria is a laggard, and Amelia is a laggard"
    else:
        return "Inconsistent assumptions"

print(check_statements())