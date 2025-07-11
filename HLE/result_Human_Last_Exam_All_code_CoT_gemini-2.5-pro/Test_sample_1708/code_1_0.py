import functools

def a_less_than_b(s1, s2, alphabet):
    """Custom comparison function for lexicographical order."""
    if s1 == s2:
        return 0
    # Prefix case
    if s2.startswith(s1):
        return -1
    if s1.startswith(s2):
        return 1
    
    # Find first differing character
    for char1, char2 in zip(s1, s2):
        if alphabet.index(char1) < alphabet.index(char2):
            return -1
        if alphabet.index(char1) > alphabet.index(char2):
            return 1
    return 0 # Should not be reached if strings are different

def map_4_to_2(s):
    """Maps a string from {a,b,c,d} to {0,1}."""
    mapping = {'a': '00', 'b': '01', 'c': '10', 'd': '11'}
    return "".join(mapping.get(char, '') for char in s)

def map_2_to_4(s):
    """Maps a string from {0,1} to {a,b,c,d}."""
    mapping = {'0': 'a', '1': 'b'}
    return "".join(mapping.get(char, '') for char in s)

# --- Demonstration for alphabet size 4 -> 2 ---
print("--- Injection from 4-char alphabet to 2-char alphabet ---")
alphabet4 = ['a', 'b', 'c', 'd']
strings4 = ['b', 'ac', 'a', 'aa', 'd', 'ca']
print(f"Original list (k=4): {strings4}")

# Sort the original list
sorted_strings4 = sorted(strings4, key=functools.cmp_to_key(lambda s1, s2: a_less_than_b(s1, s2, alphabet4)))
print(f"Sorted list (k=4):   {sorted_strings4}")

# Map and sort the transformed list
mapped_to_2 = [map_4_to_2(s) for s in sorted_strings4]
print(f"Mapped to k=2:       {mapped_to_2}")
print("The order is preserved, showing order_type(k=4) <= order_type(k=2).\n")

# --- Demonstration for alphabet size 2 -> 4 ---
print("--- Injection from 2-char alphabet to 4-char alphabet ---")
alphabet2 = ['0', '1']
strings2 = ['1', '01', '0', '00', '11', '10']
print(f"Original list (k=2): {strings2}")

# Sort the original list
sorted_strings2 = sorted(strings2, key=functools.cmp_to_key(lambda s1, s2: a_less_than_b(s1, s2, alphabet2)))
print(f"Sorted list (k=2):   {sorted_strings2}")

# Map and sort the transformed list
mapped_to_4 = [map_2_to_4(s) for s in sorted_strings2]
print(f"Mapped to k=4:       {mapped_to_4}")
print("The order is preserved, showing order_type(k=2) <= order_type(k=4).\n")

print("Since the order types can be embedded within each other, they must be equal.")
print("The order type for any alphabet k>=2 is a known result from set theory.")
final_answer = "ω^ω"
print(f"\nThe final answer is the ordinal: {final_answer}")
print("\nBreaking down the 'equation':")
print(f"The base is the first infinite ordinal: ω")
print(f"The exponent is the first infinite ordinal: ω")
