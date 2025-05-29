# Define the possible states
def is_hero(statement):
    return statement

def is_villain(statement):
    return not statement

# Chloe's statements
chloe_says_benjamin_is_hero = True
chloe_says_benjamin_is_villain = True

# Chloe is a villain, so her statements are false
chloe_is_hero = False
chloe_is_villain = True

# Benjamin's statement: "If Chloe is a hero, then Benjamin is a hero"
# This statement is only false if Chloe is a hero and Benjamin is not a hero
benjamin_statement = not (chloe_is_hero and not is_hero(True))

# Determine Benjamin's status
benjamin_is_hero = False
benjamin_is_villain = True

print(f"Benjamin is a {'hero' if benjamin_is_hero else 'villain'}, and Chloe is a {'hero' if chloe_is_hero else 'villain'}.")