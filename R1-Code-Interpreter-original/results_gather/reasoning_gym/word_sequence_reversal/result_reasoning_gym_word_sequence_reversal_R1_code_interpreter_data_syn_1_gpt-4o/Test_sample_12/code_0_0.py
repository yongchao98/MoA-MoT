# Original list of words
words = "readable, to, he, that, to, possession"

# Split the string into a list of words
words_list = words.split(", ")

# Reverse the list
reversed_list = words_list[::-1]

# Join the reversed list back into a string
reversed_string = ", ".join(reversed_list)

# Print the result
print(reversed_string)