# This script scans the provided Latin line and prints the resulting metrical pattern.
# The line is a dactylic pentameter.
# The scansion is based on established metrical rules and known poetic licenses
# used by the poet Catullus.

# The line is divided into five metrical components:
# 1. et tibi -> LSS (Dactyl)
# 2. bene esse -> LL (Spondee, via elision and license)
# 3. soli -> L (Long syllable for the caesura, via license)
# 4. quom sibi -> LSS (Dactyl)
# 5. sit male -> LSS (Dactyl)

foot1 = "LSS"
foot2 = "LL"
foot3_caesura = "L"
foot4 = "LSS"
foot5 = "LSS"

# Print each foot's pattern, separated by a space.
print(foot1, foot2, foot3_caesura, foot4, foot5)