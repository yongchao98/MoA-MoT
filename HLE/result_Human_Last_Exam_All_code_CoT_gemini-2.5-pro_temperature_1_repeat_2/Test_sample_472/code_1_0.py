# The line to be scanned is:
# "et tibi bene esse soli quom sibi sit male"
# After applying Plautine contractions (tibi -> tib', bene esse -> ben'esse, sibi -> sib'),
# the line is read as "et tib' ben'esse soli quom sib' sit male" and scanned as follows:

# Foot 1: et tib' -> L S (Trochaic substitution)
# Foot 2: ben' es-se -> S L (Iamb)
# Foot 3: -se so-li -> S L (Iamb)
# Foot 4: -li quom -> S L (Iamb)
# Foot 5: sib' sit -> S L (Iamb)
# Foot 6: ma le -> S S (Pyrrhic substitution)

foot1 = "L S"
foot2 = "S L"
foot3 = "S L"
foot4 = "S L"
foot5 = "S L"
foot6 = "S S"

# Print each foot separated by a space.
print(f"{foot1} {foot2} {foot3} {foot4} {foot5} {foot6}")
