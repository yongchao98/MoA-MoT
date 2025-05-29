# Logical check for the statements
def check_statements():
    # David's statement: "Harper is a hero and Harper is a villain."
    # This is a contradiction, so if David is a hero, this statement is false.
    # Therefore, David must be a villain.
    david_is_hero = False  # David is a villain

    # Harper's statement: "David is a hero if and only if David is a villain."
    # This is a contradiction, so if Harper is a hero, this statement is false.
    # Therefore, Harper must be a villain.
    harper_is_hero = False  # Harper is a villain

    return david_is_hero, harper_is_hero

david_is_hero, harper_is_hero = check_statements()
print(f"David is a {'hero' if david_is_hero else 'villain'}, and Harper is a {'hero' if harper_is_hero else 'villain'}.")