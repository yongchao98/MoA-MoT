# Original list of words
words = "Madrid, Smith, it"

# Split the words by comma and strip any extra spaces
word_list = [word.strip() for word in words.split(',')]

# Reverse the list of words
reversed_list = word_list[::-1]

# Join the reversed list into a string with the specified format
result = ', '.join(reversed_list)

# Print the result
print(result)