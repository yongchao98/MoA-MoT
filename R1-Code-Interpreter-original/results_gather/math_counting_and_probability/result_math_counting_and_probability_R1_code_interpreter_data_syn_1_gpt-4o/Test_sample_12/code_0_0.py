# Total number of 3-letter words
total_words = 4**3

# Number of 3-letter words without the letter A
words_without_A = 3**3

# Number of 3-letter words with at least one A
words_with_A = total_words - words_without_A

print(words_with_A)