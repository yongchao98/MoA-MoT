# The line to be scanned:
# "et tibi bene esse soli quom sibi sit male"

# This line is identified as a trochaic septenarius.
# A trochaic foot is Long-Short (L S). Common substitutions are allowed.
# The elision "bene esse" -> "ben'esse" is applied.

# Based on analysis, the line is broken down into the following metrical feet.
# Each variable represents one foot in the line.
foot_1 = "L S"    # Corresponds to "et ti"
foot_2 = "S S"    # Corresponds to "bi ben'" (from elided bene)
foot_3 = "L S"    # Corresponds to "esse"
foot_4 = "L L"    # Corresponds to "soli"
foot_5 = "L S S"  # Corresponds to "quom sibi"
foot_6 = "L S"    # Corresponds to "sit ma"
foot_7 = "S"      # Corresponds to the final syllable "-le"

# The list of feet is joined into a single string with spaces separating each foot.
final_scansion = "  ".join([foot_1, foot_2, foot_3, foot_4, foot_5, foot_6, foot_7])

# The final result is printed.
print(final_scansion)