import math
import random

def solve_language_with_magic_float():
    """
    Demonstrates how a single float can solve a complex language problem
    for a fixed input length, illustrating the principle behind the P/poly class.

    The Problem:
    - We define an arbitrary language `L` consisting of a subset of all possible
      8-bit binary strings.
    - The task is to write a function `is_in_language(s)` that returns True
      if the 8-bit string `s` is in `L`, and False otherwise.

    The "Magic Float" Solution:
    1. We create a 2^8 = 256-bit "advice string". The i-th bit is '1' if the binary
       representation of `i` is in our language `L`, and '0' otherwise.
    2. We encode this entire 256-bit advice string into a single floating-point number.
       We treat the advice string as the fractional part of a number in binary
       (e.g., 0.b_0b_1b_2...).
    3. The `is_in_language` function then uses arithmetic to extract the correct bit
       from this "magic float" to get the answer, without needing any explicit
       if/else logic or data structures for the language itself.
    """
    BIT_LENGTH = 8
    TABLE_SIZE = 2**BIT_LENGTH

    # Step 1: Define an arbitrary complex language `L` for 8-bit strings.
    # We will randomly select about half of the possible strings to be in the language.
    # The set `L` is our "non-uniform advice".
    print(f"--- Defining an arbitrary language for {BIT_LENGTH}-bit strings ---")
    language_set = set()
    for i in range(TABLE_SIZE):
        if random.random() < 0.5:
            # Format the number `i` as an 8-bit binary string (e.g., '00001101')
            binary_string = format(i, '0' + str(BIT_LENGTH) + 'b')
            language_set.add(binary_string)
    print(f"Language 'L' defined with {len(language_set)} members.\n")

    # Step 2: Encode the language's lookup table into a single "magic" float.
    # We create an "advice string" where the i-th bit is 1 if the i-th number
    # (as a binary string) is in L.
    print(f"--- Encoding the language into a single 'magic' float ---")
    magic_number = 0.0
    advice_string = ""
    for i in range(TABLE_SIZE):
        binary_string = format(i, '0' + str(BIT_LENGTH) + 'b')
        if binary_string in language_set:
            # To encode bit '1' at position i+1, we add 1 / (2^(i+1))
            magic_number += 1.0 / (2**(i + 1))
            advice_string += "1"
        else:
            advice_string += "0"

    print(f"The first 64 bits of the advice string are: {advice_string[:64]}...")
    print(f"The 'magic' number encoding this is: {magic_number:.20f}\n")
    # This number is the "weight" our hypothetical Transformer would use.


    def is_in_language(binary_string: str, magic_weight: float) -> bool:
        """
        A hypothetical function inside the Transformer.
        It uses the magic weight to decide if the string is in the language.
        """
        # Convert binary string '101' to integer 5
        index = int(binary_string, 2)

        # "Extract" the bit corresponding to this index from the magic number
        # 1. Multiply by 2^(index+1) to shift the desired bit just to the left
        #    of the decimal point.
        #    e.g., for index 2 (bit 3), 0.b_0b_1b_2... becomes b_0b_1b_2.b_3...
        shifted_value = magic_weight * (2**(index + 1))

        # 2. Taking the floor gives the integer part.
        #    e.g., floor(b_0b_1b_2.b_3...) = integer(b_0b_1b_2)
        integer_part = math.floor(shifted_value)

        # 3. If the integer part is odd, the last bit (our target bit) was 1.
        #    If it's even, the last bit was 0.
        return integer_part % 2 == 1


    # Step 3: Test the function on a few examples
    print(f"--- Testing the recognition function using the magic float ---")
    test_strings = [format(i, '0' + str(BIT_LENGTH) + 'b') for i in random.sample(range(TABLE_SIZE), 5)]

    for s in test_strings:
        is_member_actual = s in language_set
        is_member_predicted = is_in_language(s, magic_number)
        print(f"Input: {s}")
        print(f"Is '{s}' in L? (Ground Truth): {is_member_actual}")
        print(f"Is '{s}' in L? (Predicted via Magic Float): {is_member_predicted}")
        print("-" * 20)

solve_language_with_magic_float()