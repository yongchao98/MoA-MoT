# Total number of 4-letter words
total_words = 5**4

# Number of 4-letter words with no vowels
no_vowel_words = 3**4

# Number of 4-letter words with at least one vowel
at_least_one_vowel = total_words - no_vowel_words

print(at_least_one_vowel)