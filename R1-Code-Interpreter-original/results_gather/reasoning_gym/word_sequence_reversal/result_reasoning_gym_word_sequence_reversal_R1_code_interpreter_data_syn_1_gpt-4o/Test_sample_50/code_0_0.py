# Original list of words as a string
words = "allow, goals, think, Mr"

# Split the string into a list of words
words_list = words.split(", ")

# Reverse the list of words
reversed_list = list(reversed(words_list))

# Join the reversed list back into a string
reversed_string = ", ".join(reversed_list)

# Print the result
print(reversed_string)