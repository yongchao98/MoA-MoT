# 1. The numbers are derived from the syllable count of each line of the haiku.
line1_syllables = 5
line2_syllables = 7
line3_syllables = 5

# 2. The words in the haiku suggest the mathematical operations.
# "An August tempest" (5) is the starting value.
# "brings winds of change" suggests subtracting the second line's count (7).
# "A divine one yields" suggests adding the third line's count (5).
result = line1_syllables - line2_syllables + line3_syllables

# 3. The result is a number corresponding to a letter's position in the alphabet.
# We can find the letter by converting the number. A=1, B=2, C=3...
# chr(67) gives 'C', so we need to add 64 to our result.
final_letter = chr(result + 64)

# 4. Print the steps and the final answer.
print(f"The number of syllables for each line are {line1_syllables}, {line2_syllables}, and {line3_syllables}.")
print("The haiku's phrasing suggests the following equation:")
print(f"{line1_syllables} - {line2_syllables} + {line3_syllables} = {result}")
print(f"The {result}rd letter of the alphabet is the answer.")
print(f"The answer is: {final_letter}")