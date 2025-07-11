# The line to be scanned is: "et tibi bene esse soli quom sibi sit male"
# L = Long syllable
# S = Short syllable
# The meter is a Trochaic Septenarius, with feet separated by spaces.

# Foot 1: et ti (L S - Trochee)
foot_1 = "L S"

# Foot 2: bi be ne (S S S - Tribrach)
foot_2 = "S S S"

# Foot 3: es se (L S - Trochee)
foot_3 = "L S"

# Foot 4: so li (L L - Spondee)
foot_4 = "L L"

# Foot 5: quom si (L S - Trochee)
foot_5 = "L S"

# Foot 6: bi sit (S L - Iamb)
foot_6 = "S L"

# Foot 7: ma le (S L - Iamb, with final anceps syllable treated as long)
foot_7 = "S L"

# Combine the feet with spaces in between to represent the full line's scansion.
final_scansion = f"{foot_1} {foot_2} {foot_3} {foot_4} {foot_5} {foot_6} {foot_7}"

print(final_scansion)