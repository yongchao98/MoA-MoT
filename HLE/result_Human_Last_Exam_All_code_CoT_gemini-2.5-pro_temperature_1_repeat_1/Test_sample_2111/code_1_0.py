# The term to be analyzed, based on the answer to the user's question.
# The hyphen is removed for a simpler character analysis.
term = "worldeconomy"

# Define the set of vowels for easy checking.
vowels = "aeiou"

# 1. Calculate the total number of letters.
total_letters = len(term)

# 2. Count the number of vowels in the term.
vowel_count = 0
for char in term:
    if char in vowels:
        vowel_count += 1

# 3. Calculate the number of consonants.
consonant_count = total_letters - vowel_count

# 4. Print the final equation showing the relationship between the counts.
# The output will display each number in the final equation.
print(f"Analysis of the term '{term}':")
print(f"Total Letters = Vowels + Consonants")
print(f"{total_letters} = {vowel_count} + {consonant_count}")
