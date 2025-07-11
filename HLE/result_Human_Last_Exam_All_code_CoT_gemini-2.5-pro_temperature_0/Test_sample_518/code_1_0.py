# This script solves the Latin riddle by constructing the punishment word.
# The riddle uses the first syllables of four names.

# 1. "Penelopes primam..." -> The first syllable of Penelope.
syllable1 = "Penelope"[:2].upper()

# 2. "...Didonis prima..." -> The first syllable of Dido.
syllable2 = "Dido"[:2].upper()

# 3. "...primam Cadmi..." -> The first syllable of Cadmus.
syllable3 = "Cadmus"[:2].upper()

# 4. "...syllaba prima Remi" -> The first syllable of Remus.
syllable4 = "Remus"[:2].upper()

# The riddle combines these syllables to form the punishment word.
# "quodque fit ex illis..." -> "and what is made from those things..."
punishment_word = syllable1 + syllable2 + syllable3 + syllable4

print("The punishment is a word constructed from four parts derived from the riddle:")
print(f"Part 1, from 'Penelopes': {syllable1}")
print(f"Part 2, from 'Didonis': {syllable2}")
print(f"Part 3, from 'Cadmi': {syllable3}")
print(f"Part 4, from 'Remi': {syllable4}")
print("\nCombining these parts forms the final word:")
print(f"{syllable1} + {syllable2} + {syllable3} + {syllable4} = {punishment_word.lower()}")
print(f"\nThe punishment for stealing from the garden is the Latin verb '{punishment_word.lower()}'.")