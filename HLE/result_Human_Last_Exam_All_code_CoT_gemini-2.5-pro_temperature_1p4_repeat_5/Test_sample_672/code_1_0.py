# This script performs scansion on the given line of poetry.
# The line is: "And in the letter, my cousin mentions a piece of advice"
# "/" represents a stressed syllable.
# "x" represents an unstressed syllable.

# Step 1: Analyze the syllables and their natural stress.
# "And": unstressed (x)
# "in": unstressed (x)
# "the": unstressed (x)
# "let-ter": stressed, unstressed (/x)
# "my": unstressed (x)
# "cou-sin": stressed, unstressed (/x)
# "men-tions": stressed, unstressed (/x)
# "a": unstressed (x)
# "piece": stressed (/)
# "of": unstressed (x)
# "ad-vice": unstressed, stressed (x/)

# Step 2: Combine the stress markers for all syllables.
scansion_result = "xxx/xx/x/xx/xx/"

# Step 3: Print the final scansion string.
print(scansion_result)