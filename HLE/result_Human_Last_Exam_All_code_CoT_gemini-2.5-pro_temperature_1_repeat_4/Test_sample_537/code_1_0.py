def calculate_value(word):
  """Calculates the value of a word based on the number of vowels and consonants."""
  vowels = "aeiou"
  vowel_count = 0
  for char in word:
    if char in vowels:
      vowel_count += 1
  consonant_count = len(word) - vowel_count
  return vowel_count * consonant_count

# The list of words from the puzzle
words = ["utpshtheas", "fkuhu", "numhkatnatu", "anhaeketn"]

# Calculate the values for the first three words
val1 = calculate_value(words[0])
val2 = calculate_value(words[1])
val3 = calculate_value(words[2])

# Print the equation that shows the pattern
# The instruction is to "output each number in the final equation!"
print(val1, "+", val2, "+", 1, "=", val3)

# Calculate and print the value for the last word
final_value = calculate_value(words[3])
print(final_value)