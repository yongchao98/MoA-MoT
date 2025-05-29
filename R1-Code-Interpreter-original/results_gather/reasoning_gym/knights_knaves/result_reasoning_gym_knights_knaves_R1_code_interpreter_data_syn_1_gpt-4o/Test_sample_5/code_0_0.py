# Define the logic for the statements
def check_statements():
    # Assume Henry is a hero
    henry_is_hero = True
    ethan_is_hero = henry_is_hero  # According to Henry's statement
    ethan_statement = not henry_is_hero  # Ethan's statement

    # Check if the assumptions lead to a consistent scenario
    if ethan_is_hero == ethan_statement:
        return "Henry is a hero, and Ethan is a villain"
    else:
        return "Henry is a villain, and Ethan is a hero"

print(check_statements())