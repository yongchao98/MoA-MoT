def check_statements():
    # Assume Amelia is a hero (truth-teller)
    amelia_is_hero = True
    jack_is_hero = False  # Since Amelia's statement implies Jack is a villain

    # Check Amelia's statement: If Jack is a hero, then Jack is a villain
    amelia_statement = not jack_is_hero or not jack_is_hero  # This should be True if Amelia is a hero

    # Check Jack's statement: Amelia is a villain if and only if Amelia is a hero
    jack_statement = (not amelia_is_hero) == amelia_is_hero  # This should be False if Jack is a villain

    return amelia_statement, not jack_statement

amelia_statement, jack_statement = check_statements()
print(amelia_statement, jack_statement)