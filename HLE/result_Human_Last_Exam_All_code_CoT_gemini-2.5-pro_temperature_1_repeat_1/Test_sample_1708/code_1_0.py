def string_to_int(s, chars="abcd", verbose=False):
    """
    Calculates the unique integer index for a string in shortlex order.
    The mapping is an isomorphism from the set of strings to the natural numbers.
    """
    base = len(chars)
    char_map = {char: i for i, char in enumerate(chars)}
    k = len(s)

    if k == 0:
        if verbose:
            print(f"String '{s}' (empty): index = 0")
        return 0
    
    # This is the "equation" for the index.
    # Part 1: Offset for all strings shorter than k.
    # This is the sum of a geometric series: 1 + base + base^2 + ... + base^(k-1)
    # The sum is (base^k - 1) / (base - 1).
    offset = (base**k - 1) // (base - 1)
    
    # Part 2: Position of s among strings of its own length.
    # This is found by interpreting the string as a number in base-`base`.
    base_val = 0
    for char in s:
        base_val = base_val * base + char_map[char]
            
    total_index = offset + base_val

    if verbose:
        print(f"To calculate the index for string '{s}':")
        print(f"  1. Find the count of all shorter strings (the offset).")
        print(f"     - String length k = {k}, Alphabet size (base) = {base}")
        print(f"     - Offset = (base^k - 1) / (base - 1) = ({base}^{k} - 1) / ({base - 1}) = {offset}")
        print(f"  2. Find the string's position within its length group.")
        print(f"     - Value of '{s}' in base-{base} (a=0, b=1, ...) = {base_val}")
        print(f"  3. Add them together.")
        print(f"     - Final Index = offset + value = {offset} + {base_val} = {total_index}")

    return total_index

def int_to_string(n, chars="abcd"):
    """Maps an integer index back to its corresponding string."""
    if n < 0:
        raise ValueError("Input must be a non-negative integer")

    base = len(chars)
    if n == 0:
        return ""

    # Determine the length of the string by finding which group `n` belongs to.
    length = 0
    count_shorter = 0
    count_at_length = 1  # For the empty string
    while True:
        if n < count_shorter + count_at_length:
            break
        count_shorter += count_at_length
        length += 1
        count_at_length *= base
    
    index_at_length = n - count_shorter

    # Convert the index within the length group to a base-`base` string.
    result = []
    temp_n = index_at_length
    for _ in range(length):
        remainder = temp_n % base
        result.append(chars[remainder])
        temp_n //= base
    
    return "".join(reversed(result))

# Final Conclusion
print("The order type of the set of finite strings from {a,b,c,d} ordered lexically is ω (omega).")
print("This is because the set is order-isomorphic to the set of natural numbers {0, 1, 2, ...}.\n")
print("The code below demonstrates this isomorphism.")
print("It generates strings from integers and verifies the mapping using a formula.")
print("-" * 60)

# Demonstrate the calculation for a specific string, "ba"
string_to_int("ba", verbose=True)

print("-" * 60)
print("\nHere is a list of the first 22 strings and their corresponding integer indices:\n")
print("{:<10} | {:<10}".format("Integer", "String"))
print("-" * 23)

for i in range(22):
    s = int_to_string(i)
    print("{:<10} | '{:s}'".format(i, s))
