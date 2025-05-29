def check_statements():
    # Case 1: Zoey is a sage
    zoey_is_sage = True
    riley_is_fool = True
    riley_statement_case1 = (zoey_is_sage or not riley_is_fool)  # Riley's statement should be false
    case1_valid = riley_statement_case1 == False

    # Case 2: Zoey is a fool
    zoey_is_fool = True
    riley_is_sage = True
    riley_statement_case2 = (not zoey_is_fool or riley_is_sage)  # Riley's statement should be true
    case2_valid = riley_statement_case2 == True

    if case1_valid:
        return "Zoey is a sage, and Riley is a fool"
    elif case2_valid:
        return "Zoey is a fool, and Riley is a sage"
    else:
        return "No consistent solution"

print(check_statements())