# The plan is to first identify all words in the sentence.
# Then, for each word, determine the number of syllables and the position of the stress.
# A simple way to count syllables is to count the vowels (а, е, ё, и, о, у, ы, э, ю, я).
# Words with only one vowel (one syllable) will be excluded.
# For the remaining words, we identify the stress position. This requires linguistic knowledge of Russian.
# We will filter for words where the stress is not on the first syllable.
# Finally, we will print the resulting words, separated by commas, in their order of appearance.

# 1. Text: "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."
# 2. Word analysis:
#    - Шашлык: 2 syllables (Ша-шлы́к). Stress on syllable 2. -> Keep.
#    - запах: 2 syllables (за́-пах). Stress on syllable 1. -> Discard.
#    - горелым: 3 syllables (го-ре́-лым). Stress on syllable 2. -> Keep.
#    - и: 1 syllable. -> Discard.
#    - вскоре: 2 syllables (вско́-ре). Stress on syllable 1. -> Discard.
#    - прибежал: 3 syllables (при-бе-жа́л). Stress on syllable 3. -> Keep.
#    - сосед: 2 syllables (со-се́д). Stress on syllable 2. -> Keep.
#    - видно: 2 syllables (ви́д-но). Stress on syllable 1. -> Discard.
#    - запах: 2 syllables (за́-пах). Stress on syllable 1. -> Discard.
#    - дошёл: 2 syllables (до-шёл). Stress on syllable 2. -> Keep.
#    - и: 1 syllable. -> Discard.
#    - до: 1 syllable. -> Discard.
#    - его: 2 syllables (е-го́). Stress on syllable 2. -> Keep.
#    - дома: 2 syllables (до́-ма). In this context, stress is on syllable 1. -> Discard.

# 3. Final list of words based on the analysis:
words_with_non_first_syllable_accent = ['Шашлык', 'горелым', 'прибежал', 'сосед', 'дошёл', 'его']

# 4. Print the result as a comma-separated string.
print(", ".join(words_with_non_first_syllable_accent))