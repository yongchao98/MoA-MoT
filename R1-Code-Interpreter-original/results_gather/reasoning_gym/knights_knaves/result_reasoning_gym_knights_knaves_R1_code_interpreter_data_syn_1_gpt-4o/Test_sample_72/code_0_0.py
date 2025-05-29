def check_statements():
    # Scenario 1: Ethan is a hero
    ethan_hero = True
    liam_hero = ethan_hero  # Ethan's statement: Ethan is a hero if and only if Liam is a hero
    if liam_hero:
        ethan_hero = True  # Liam's statement: If Liam is a hero then Ethan is a hero
    scenario_1_consistent = ethan_hero and liam_hero

    # Scenario 2: Ethan is a villain
    ethan_hero = False
    liam_hero = not ethan_hero  # Ethan's statement: Ethan is a hero if and only if Liam is a hero
    if liam_hero:
        ethan_hero = True  # Liam's statement: If Liam is a hero then Ethan is a hero
    scenario_2_consistent = not ethan_hero and not liam_hero

    return scenario_1_consistent, scenario_2_consistent

print(check_statements())