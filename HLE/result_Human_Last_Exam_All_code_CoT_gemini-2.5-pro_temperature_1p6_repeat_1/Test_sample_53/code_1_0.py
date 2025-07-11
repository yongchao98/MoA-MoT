import sys

# This script simulates the phonological development of a PIE root into Middle English.
# The transformations are based on established sound laws in historical linguistics.

# Step 1: From Proto-Indo-European to Proto-Germanic
pie_root = "*kʷeys-"
pie_causative_stem = "*kʷoyéye-"
pgmc_infinitive = "*gʷajjaną"

print("This script will trace the hypothetical inheritance of a PIE root into Middle English.")
print("---")
print(f"1. Start with Proto-Indo-European root: {pie_root} ('to see, to heed').")
print(f"2. Form the o-grade causative stem, meaning 'to cause to see' or 'to show': {pie_causative_stem}.")
print(f"3. Apply sound changes to get to Proto-Germanic '{pgmc_infinitive}':")
print("   - Grimm's Law changes *kʷ > *hʷ.")
print("   - Verner's Law then voices *hʷ > *gʷ because the root syllable was unstressed.")
print("   - PIE *o becomes PGmc *a.")
print("   - The causative suffix *-éye- becomes *-ja- in a Class 1 weak verb.")

# Step 2: From Proto-Germanic to Old English
oe_infinitive = "cweċċan"

print("\n---")
print(f"4. Trace from Proto-Germanic '{pgmc_infinitive}' to the Old English infinitive '{oe_infinitive}':")
print("   - I-umlaut: The 'j' in *-janą fronts the root vowel 'a' to 'e'. (*gʷajjaną > *gʷejjaną)")
print("   - Consonant Shift: PGmc *gʷ becomes OE 'cw' before a front vowel like 'e'. (*gʷejjaną > cwejjan...)")
print("   - Palatalization: The '-jj-' sound becomes the affricate sound '-ċċ-' (/tʃtʃ/). (cwejjan > cweċċan)")

# Step 3: From Old English to Middle English
oe_3s_present = "cweċċeþ"
me_3s_present = "quecheth"

print("\n---")
print(f"5. Form the 3rd person singular present in Old English: add the suffix '-eþ'.")
print(f"   - {oe_infinitive} -> {oe_3s_present} ('he shows')")
print(f"6. Trace from Old English '{oe_3s_present}' to the Middle English form '{me_3s_present}':")
print("   - Spelling: The initial 'cw' is regularly spelled 'qu' in Middle English.")
print("   - Simplification: The double affricate 'ċċ' /tʃtʃ/ simplifies to 'ch' /tʃ/.")
print("   - Ending: The verbal ending '-eþ' becomes '-eth'.")

# Final equation and answer
print("\n---")
print("Final Transformation:")
print(f"Old English '{oe_3s_present}' -> Middle English '{me_3s_present}'")
print(f"\nThe predicted third person singular present verbal form in Middle English is:")
print(me_3s_present)

# Set the final answer for the system.
# We use a file object that the system can read.
try:
    with open('output.txt', 'w') as f:
        f.write(f"<<<{me_3s_present}>>>")
except Exception:
    # If file writing fails, print to stdout for compatibility
    sys.stdout.write(f"\n<<<{me_3s_present}>>>")
