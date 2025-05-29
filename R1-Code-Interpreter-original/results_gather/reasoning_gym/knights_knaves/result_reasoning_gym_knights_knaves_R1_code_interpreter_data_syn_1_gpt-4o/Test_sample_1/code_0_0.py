def check_statements():
    # Scenario 1: Elizabeth is a saint
    elizabeth_is_saint = True
    ethan_is_sinner = True
    ethan_statement_true = not ethan_is_sinner or not elizabeth_is_saint
    if elizabeth_is_saint and ethan_is_sinner and not ethan_statement_true:
        return "Elizabeth is a saint, and Ethan is a sinner"

    # Scenario 2: Elizabeth is a sinner
    elizabeth_is_saint = False
    ethan_is_sinner = False
    ethan_statement_true = not ethan_is_sinner or not elizabeth_is_saint
    if not elizabeth_is_saint and not ethan_is_sinner and ethan_statement_true:
        return "Elizabeth is a sinner, and Ethan is a saint"

    return "No consistent scenario found"

print(check_statements())