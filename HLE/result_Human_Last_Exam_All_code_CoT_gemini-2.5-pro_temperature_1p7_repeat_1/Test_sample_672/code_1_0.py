# This script performs scansion on the given line of poetry.
# Scansion is the process of marking the stressed and unstressed syllables.
# '/' is used for a stressed syllable.
# 'x' is used for an unstressed syllable.

# The line of poetry is:
# "And in the letter, my cousin mentions a piece of advice"

# The analysis proceeds as follows:
# "And" - unstressed (x)
# "in" - unstressed (x)
# "the" - unstressed (x)
# "let-ter," - stressed, unstressed (/x)
# "my" - unstressed (x)
# "cou-sin" - stressed, unstressed (/x)
# "men-tions" - stressed, unstressed (/x)
# "a" - unstressed (x)
# "piece" - stressed (/)
# "of" - unstressed (x)
# "ad-vice" - unstressed, stressed (x/)
scansion = "x" + "x" + "x" + "/" + "x" + "x" + "/" + "x" + "/" + "x" + "x" + "/" + "x" + "x" + "/"

# Print the final scansion string
print(scansion)