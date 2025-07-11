# This script prints the scansion of the Latin line:
# "et tibi bene esse soli quom sibi sit male"
# L = Long syllable, S = Short syllable.

# My analysis shows the meter is Iambic Septenarius, and the scansion is as follows:
foot1 = "L S"  # et ti (Trochee substitution)
foot2 = "S L"  # bi ben (Iamb)
foot3 = "L S"  # es se (Trochee substitution)
foot4 = "L L"  # so li (Spondee substitution)
foot5 = "L S"  # quom si (Trochee substitution)
foot6 = "S L"  # bi sit (Iamb)
foot7 = "S S"  # ma le (Pyrrhic substitution)

# Joining the feet with a space in between to form the complete line's scansion.
# The instruction to "output each number" is interpreted as providing the
# final sequence of L and S characters that make up the scansion.
full_scansion = " ".join([foot1, foot2, foot3, foot4, foot5, foot6, foot7])

# Print the final result.
print(full_scansion)