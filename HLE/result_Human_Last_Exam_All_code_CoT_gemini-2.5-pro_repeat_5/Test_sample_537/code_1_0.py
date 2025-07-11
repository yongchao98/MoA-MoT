def count_unique_vowels(word):
    """Counts the number of unique vowels in a word."""
    vowels = "aeiou"
    unique_vowels_found = set()
    for char in word:
        if char in vowels:
            unique_vowels_found.add(char)
    return len(unique_vowels_found)

# The given words
word1 = "utpshtheas"
word2 = "fkuhu"
word3 = "numhkatnatu"
word4 = "anhaeketn"

# Step 1: Decode each word by counting its unique vowels
n1 = count_unique_vowels(word1)
n2 = count_unique_vowels(word2)
n3 = count_unique_vowels(word3)
n4 = count_unique_vowels(word4)

# Step 2: The pattern is n1 - n2 = n3. Let's verify.
# 3 - 1 = 2, which is correct.

# Step 3: Apply the pattern to the next numbers in the sequence (n3 and n4)
# to find the result for the last word.
result = n3 - n4

# Step 4: Print the final equation with each number.
print(f"{n3} - {n4} = {result}")

# The decoded value of the last word is the result of the equation.
# final_answer = result