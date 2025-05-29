def is_sarcastic(statement):
    negative_consequences = "no water and no electricity"
    if negative_consequences in statement and "smart move" in statement:
        return True
    return False

options = {
    "A": "Without the river there's no water and no electricity. Blocking it is such a smart move!",
    "B": "Without the river there's no water and no electricity. Blocking it is such a terrible move!"
}

sarcastic_option = None
for option, statement in options.items():
    if is_sarcastic(statement):
        sarcastic_option = option
        break

print(sarcastic_option)