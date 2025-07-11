# This script performs scansion on the line of poetry.
# The line is: "And in the letter, my cousin mentions a piece of advice"
# The scansion uses "/" for stressed syllables and "x" for unstressed syllables.

# Step 1: Syllabification and stress analysis
# "And" -> x (unstressed)
# "in" -> x (unstressed)
# "the" -> x (unstressed)
# "let-ter," -> /x (stressed, unstressed)
# "my" -> x (unstressed)
# "cou-sin" -> /x (stressed, unstressed)
# "men-tions" -> /x (stressed, unstressed)
# "a" -> x (unstressed)
# "piece" -> / (stressed)
# "of" -> x (unstressed)
# "ad-vice" -> x/ (unstressed, stressed)

# Step 2: Combine the stress markers into the final scansion string.
scansion_result = "xxx/xx/x/xx/xx/"

# Step 3: Print the final result without any extra text or formatting.
print(scansion_result)