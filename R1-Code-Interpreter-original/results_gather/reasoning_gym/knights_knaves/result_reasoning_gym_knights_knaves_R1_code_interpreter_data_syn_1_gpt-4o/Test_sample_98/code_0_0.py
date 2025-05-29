def check_statements():
    # Scenario 1: Sebastian is an altruist
    sebastian_altruist = True
    ethan_egoist = True
    if sebastian_altruist and ethan_egoist:
        sebastian_statement = True
        ethan_statement = False
    else:
        sebastian_statement = False
        ethan_statement = True

    # Check if both statements are consistent
    scenario_1_valid = sebastian_statement and not ethan_statement

    # Scenario 2: Sebastian is an egoist
    sebastian_altruist = False
    ethan_egoist = False
    if not sebastian_altruist and not ethan_egoist:
        sebastian_statement = False
        ethan_statement = True
    else:
        sebastian_statement = True
        ethan_statement = False

    # Check if both statements are consistent
    scenario_2_valid = not sebastian_statement and ethan_statement

    return scenario_1_valid, scenario_2_valid

scenario_1_valid, scenario_2_valid = check_statements()
print(scenario_1_valid, scenario_2_valid)