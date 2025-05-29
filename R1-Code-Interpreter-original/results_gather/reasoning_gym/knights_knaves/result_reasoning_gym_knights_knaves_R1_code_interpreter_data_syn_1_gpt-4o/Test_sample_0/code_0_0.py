def check_statements():
    # Scenario 1: Scarlett is a saint
    scarlett_saint = True
    victoria_sinner = True  # Based on Scarlett's statement
    victoria_statement = victoria_sinner or not scarlett_saint
    if scarlett_saint and not victoria_statement:
        return "Contradiction in Scenario 1"

    # Scenario 2: Scarlett is a sinner
    scarlett_saint = False
    victoria_saint = True  # Based on Scarlett's statement being false
    victoria_statement = victoria_saint or not scarlett_saint
    if not scarlett_saint and victoria_statement:
        return "Scarlett is a sinner, and Victoria is a saint"

    return "No valid scenario found"

print(check_statements())