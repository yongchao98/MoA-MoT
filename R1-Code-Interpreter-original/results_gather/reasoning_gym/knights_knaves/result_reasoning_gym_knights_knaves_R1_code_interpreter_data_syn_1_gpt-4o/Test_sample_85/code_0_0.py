def check_statements():
    # Possible states: True for altruist, False for egoist
    for ava in [True, False]:
        for victoria in [True, False]:
            # Ava's statement: "Victoria is an altruist if and only if Ava is an altruist."
            ava_statement = (victoria == ava)
            # Victoria's statement: "Ava is an egoist."
            victoria_statement = not ava
            
            # Check if the statements are consistent with their nature
            if (ava == ava_statement) and (victoria == victoria_statement):
                return ava, victoria

ava_is_altruist, victoria_is_altruist = check_statements()
print(f"Ava is an {'altruist' if ava_is_altruist else 'egoist'}, and Victoria is an {'altruist' if victoria_is_altruist else 'egoist'}.")