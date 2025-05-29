import string

# Known incorrect letters
incorrect_letters = set("QSYMBURHWPCAGOJ")

# All possible letters
all_letters = set(string.ascii_uppercase)

# Possible letters are those not in the incorrect list
possible_letters = all_letters - incorrect_letters

# We know 'I' is one of the letters, but in the wrong position in guess 8
# So, the other letter must be from the possible letters
for letter in possible_letters:
    if letter != 'I':
        print(["3", "0", "I", letter])
        break