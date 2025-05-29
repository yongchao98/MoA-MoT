def check_statements():
    # Case 1: Assume Joseph is a saint
    joseph_is_saint = True
    lily_is_sinner = True  # From Joseph's statement if he is a saint
    lily_statement = not joseph_is_saint  # Lily says Joseph is a sinner

    if joseph_is_saint and lily_is_sinner and not lily_statement:
        return "Joseph is a saint, and Lily is a sinner"

    # Case 2: Assume Joseph is a sinner
    joseph_is_saint = False
    lily_is_saint = True  # From Lily's statement if she is a saint
    joseph_statement = not (joseph_is_saint and not lily_is_saint)

    if not joseph_is_saint and lily_is_saint and joseph_statement:
        return "Joseph is a sinner, and Lily is a saint"

    return "No consistent solution found"

print(check_statements())