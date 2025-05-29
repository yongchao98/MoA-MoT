def check_statements():
    # Scenario 1: Emma is an angel
    emma_is_angel = True
    jackson_is_angel = emma_is_angel  # Based on Emma's statement
    jackson_statement = emma_is_angel  # Jackson's statement

    scenario_1_valid = jackson_is_angel == jackson_statement

    # Scenario 2: Emma is a devil
    emma_is_angel = False
    jackson_is_angel = emma_is_angel  # Based on Emma's statement
    jackson_statement = emma_is_angel  # Jackson's statement

    scenario_2_valid = jackson_is_angel != jackson_statement

    return scenario_1_valid, scenario_2_valid

print(check_statements())