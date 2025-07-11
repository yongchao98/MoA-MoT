# This script performs scansion on the line of poetry.
# The process involves analyzing the stress pattern of each syllable in the line.
# "/" is used for a stressed syllable.
# "x" is used for an unstressed syllable.

# The line is: "And in the letter, my cousin mentions a piece of advice"

# Breakdown of the scansion:
# "And"      -> x (unstressed)
# "in"       -> x (unstressed)
# "the"      -> x (unstressed)
# "let-ter"  -> /x (stressed, unstressed)
# "my"       -> x (unstressed)
# "cou-sin"  -> /x (stressed, unstressed)
# "men-tions"-> /x (stressed, unstressed)
# "a"        -> x (unstressed)
# "piece"    -> / (stressed)
# "of"       -> x (unstressed)
# "ad-vice"  -> x/ (unstressed, stressed)

# Combining these marks gives the final scansion string.
final_scansion = "xxx/xx/x/xx/x/"

print(final_scansion)