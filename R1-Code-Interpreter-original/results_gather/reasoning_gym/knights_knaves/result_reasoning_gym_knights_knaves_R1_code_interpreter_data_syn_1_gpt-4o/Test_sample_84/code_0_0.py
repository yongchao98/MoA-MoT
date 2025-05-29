def check_statements():
    # Assume Victoria is a sage
    victoria_is_sage = True
    grace_is_fool = True  # Based on Victoria's statement if she is a sage

    # Check Grace's statement
    grace_statement = victoria_is_sage and not grace_is_fool  # Grace says both are sages

    # If Grace's statement is false, Grace is a fool
    if not grace_statement:
        print("Victoria is a sage, and Grace is a fool")
    else:
        print("Victoria is a fool, and Grace is a sage")

check_statements()