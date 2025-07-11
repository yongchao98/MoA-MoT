# The clue "Кома" is a cryptogram based on chemical elements.
# The letters К, О, М, А are the initials of elements.
# We need to find their atomic numbers and sum them.

# Atomic number for Калий (Potassium, K)
potassium = 19
# Atomic number for Oxygen (O)
oxygen = 8
# Atomic number for Магний (Magnesium, Mg)
magnesium = 12
# Atomic number for Аргон (Argon, Ar)
argon = 18

# The sum of the atomic numbers reveals the administrative code of the location.
location_code = potassium + oxygen + magnesium + argon

print("The clue 'Кома' can be decoded by summing the atomic numbers of the elements its letters represent:")
print(f"К (Potassium): {potassium}")
print(f"О (Oxygen): {oxygen}")
print(f"М (Magnesium): {magnesium}")
print(f"А (Argon): {argon}")
print("\nThe final equation is:")
print(f"{potassium} + {oxygen} + {magnesium} + {argon} = {location_code}")
print(f"\nThe sum, {location_code}, is the OKATO administrative code for Perm Krai.")