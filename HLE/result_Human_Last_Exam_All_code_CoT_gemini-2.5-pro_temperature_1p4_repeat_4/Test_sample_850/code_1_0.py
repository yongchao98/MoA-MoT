import textwrap

def print_wrapped(text):
    """Helper function to print text with wrapping for better readability."""
    print(textwrap.fill(text, width=80))

# --- Define Crystal Class Sets ---

# All 32 point groups
all_crystal_classes = {
    "1", "-1", "2", "m", "2/m", "222", "mm2", "mmm", "4", "-4", "4/m",
    "422", "4mm", "-42m", "4/mmm", "3", "-3", "32", "3m", "-3m", "6",
    "-6", "6/m", "622", "6mm", "-6m2", "6/mmm", "23", "m-3", "432",
    "-43m", "m-3m"
}

# 1. Optically Active Classes: These are the 11 CHIRAL classes, which lack mirror planes and inversion centers.
optically_active_classes = {"1", "2", "3", "4", "6", "222", "32", "422", "622", "23", "432"}

# 2. Achiral Classes: These are all classes that are NOT chiral.
achiral_classes = all_crystal_classes - optically_active_classes

# 3. Non-polar Classes: These are the 11 centrosymmetric classes, which possess a center of inversion.
# All non-polar classes are inherently achiral.
non_polar_classes = {"-1", "2/m", "mmm", "4/m", "-3", "6/m", "4/mmm", "-3m", "6/mmm", "m-3", "m-3m"}

# --- Analysis ---

print("--- Analysis of Crystal Properties for Optical Activity ---")
print_wrapped("\nThe question asks for crystal classes that are simultaneously achiral, non-polar, and have the correct symmetry for optical activity.")
print("-" * 50)

print_wrapped("\nStep 1: The primary requirement for optical activity is that the crystal must be CHIRAL. Chiral crystal classes are those without a mirror plane or a center of inversion.")
print(f"Optically Active (Chiral) Classes ({len(optically_active_classes)} total):")
print(f"-> {sorted(list(optically_active_classes))}")

print_wrapped("\nStep 2: The question specifies that the classes must be ACHIRAL. An achiral class, by definition, MUST contain a mirror plane or a center of inversion. This is the direct opposite of the requirement for optical activity.")
print(f"Achiral Classes ({len(achiral_classes)} total):")
print(f"-> {sorted(list(achiral_classes))}")

print_wrapped("\nStep 3: The question also specifies the classes must be NON-POLAR. Non-polar classes are those that have a center of inversion. Note that all non-polar classes are a subset of the achiral classes.")
print(f"Non-Polar (Centrosymmetric) Classes ({len(non_polar_classes)} total):")
print(f"-> {sorted(list(non_polar_classes))}")
print("-" * 50)


# --- Conclusion ---
# Find the intersection of the three sets as requested.
# Logically, the intersection between optically_active_classes and achiral_classes is empty.
result = achiral_classes.intersection(non_polar_classes).intersection(optically_active_classes)

print("\n--- Final Result ---")
print_wrapped("We are looking for the 'equation' or intersection: (Achiral) AND (Non-Polar) AND (Optically Active).")
print("\nFinal 'Equation':")
print(f"Set of Achiral Classes ∩ Set of Non-Polar Classes ∩ Set of Optically Active Classes = ?")
print("\nResult:")
if not result:
    print_wrapped("The resulting set is empty. There are no crystal classes that satisfy these contradictory conditions. A crystal cannot be both achiral and optically active.")
    print("Final equation answer: {}")
else:
    # This block will not be reached
    print(f"Resulting classes: {result}")
