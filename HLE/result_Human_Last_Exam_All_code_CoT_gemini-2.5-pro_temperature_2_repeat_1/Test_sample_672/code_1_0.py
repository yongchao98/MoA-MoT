# The task is to perform scansion on a line of poetry.
# The method is to identify the stress pattern syllable by syllable.
# Stressed syllables are marked with '/' and unstressed with 'x'.

# The line is: "And in the letter, my cousin mentions a piece of advice"

# Step 1: Analyze the stress pattern of each word or phrase.
# "And in the" -> Unstressed function words. Pattern: "xxx"
# "letter" -> Noun with stress on the first syllable. Pattern: "/x"
# "my" -> Unstressed pronoun. Pattern: "x"
# "cousin" -> Noun with stress on the first syllable. Pattern: "/x"
# "mentions" -> Verb with stress on the first syllable. Pattern: "/x"
# "a" -> Unstressed article. Pattern: "x"
# "piece" -> Stressed noun. Pattern: "/"
# "of" -> Unstressed preposition. Pattern: "x"
# "advice" -> Noun with stress on the second syllable. Pattern: "x/"

# Step 2: Combine the patterns in sequence to form the full scansion.
full_scansion = "xxx" + "/x" + "x" + "/x" + "/x" + "x" + "/" + "x" + "x/"

# Step 3: Print the final result without any spaces.
print(full_scansion)