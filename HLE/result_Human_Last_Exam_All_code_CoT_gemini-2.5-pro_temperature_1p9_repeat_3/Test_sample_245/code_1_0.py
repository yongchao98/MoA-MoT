# Plan:
# 1. The original Russian sentence is: "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."
# 2. We identify words with more than one syllable and where the accent is not on the first syllable.
#
#    - "Шашлы́к": 2 syllables, accent on the 2nd. -> Keep.
#    - "за́пах": 2 syllables, accent on the 1st. -> Exclude.
#    - "горе́лым": 3 syllables, accent on the 2nd. -> Keep.
#    - "и": 1 syllable. -> Exclude.
#    - "вско́ре": 2 syllables, accent on the 1st. -> Exclude.
#    - "прибежа́л": 3 syllables, accent on the 3rd. -> Keep.
#    - "сосе́д": 2 syllables, accent on the 2nd. -> Keep.
#    - "ви́дно": 2 syllables, accent on the 1st. -> Exclude.
#    - "за́пах": (already analyzed) -> Exclude.
#    - "дошёл": 2 syllables, accent on the 2nd (the letter 'ё' is always stressed). -> Keep.
#    - "и": (already analyzed) -> Exclude.
#    - "до": 1 syllable. -> Exclude.
#    - "его́": 2 syllables, accent on the 2nd. -> Keep.
#    - "до́ма": 2 syllables, accent on the 1st. -> Exclude.
#
# 3. The final list of words in order of appearance is: ["Шашлык", "горелым", "прибежал", "сосед", "дошёл", "его"].
# 4. This script will join these words with a comma and print the result.

selected_words = ["Шашлык", "горелым", "прибежал", "сосед", "дошёл", "его"]
output_string = ",".join(selected_words)
print(output_string)