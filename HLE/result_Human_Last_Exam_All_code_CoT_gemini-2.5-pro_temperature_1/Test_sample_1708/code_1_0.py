from fractions import Fraction

def string_to_fraction(s, alphabet):
    """
    Converts a string to a fraction based on its lexicographical position.
    This is an order-isomorphism. It maps strings over an alphabet of size k
    to base-(k+1) fractions using digits {1, ..., k}.
    """
    # The empty string is the first element, mapped to 0.
    if not s:
        return Fraction(0)

    # The alphabet has k characters. We use a base-(k+1) representation.
    k = len(alphabet)
    base = k + 1
    
    # Map characters 'a', 'b', 'c', 'd' to digits 1, 2, 3, 4.
    char_map = {char: i + 1 for i, char in enumerate(alphabet)}

    num = Fraction(0)
    power_of_base = Fraction(1)

    for char in s:
        digit = char_map.get(char)
        if digit is None:
            raise ValueError(f"String contains character not in alphabet: {char}")
        power_of_base /= base
        num += digit * power_of_base

    return num

# --- Main Execution ---

# The set of characters for the strings.
alphabet = ['a', 'b', 'c', 'd']

# A sample of strings to demonstrate the ordering.
strings = [
    "", "a", "aa", "ab", "b", "c", "d", "ad", "ac", "ba", "ddd", "dd"
]

# Sort the strings lexicographically to see the correct order.
sorted_strings = sorted(list(set(strings)))

print("Demonstration of the order-preserving map from strings to fractions:")
print("-" * 65)
print(f"{'Lexicographically Ordered String':<35} | {'Mapped Fraction':<15} | {'Decimal'}")
print("-" * 65)

for s_val in sorted_strings:
    # This mapping preserves the lexicographical order.
    frac = string_to_fraction(s_val, alphabet)
    print(f"{s_val!r:<35} | {str(frac):<15} | {float(frac):.10f}")

print("-" * 65)
print("\nThe order type of this set, which we can call tau (τ), is the unique")
print("non-empty order type that satisfies the following recursive equation:")

# The equation that defines the order type is tau = 1 + 4 * tau.
# The numbers in this equation are 1 and 4.
num1 = 1
num2 = len(alphabet)

print(f"\nτ = {num1} + {num2} ⋅ τ")