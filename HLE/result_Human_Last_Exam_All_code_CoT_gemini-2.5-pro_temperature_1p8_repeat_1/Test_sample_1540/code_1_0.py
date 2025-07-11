# This script calculates the memory usage for the described C program.

# 1. Memory for the first string.
# The string can be up to 100 characters long. In C, we need one extra
# byte for the null terminator '\0'. The size of a char is 1 byte.
mem_string_1 = 101 * 1

# 2. Memory for the loop index.
# An integer variable is needed to keep track of the position in the string.
# On a typical system for competitive programming, an int takes 4 bytes.
mem_index_variable = 4

# 3. Memory for the character from the second string.
# To read the second string character by character without storing it, we can use
# getchar(). getchar() returns an int, which takes 4 bytes.
mem_char_variable = 4

# 4. Total memory calculation.
# The total memory is the sum of the memory for these three variables.
m = mem_string_1 + mem_index_variable + mem_char_variable

# Output the breakdown of the calculation as requested.
print(f"Memory for the first string (char s1[101]): {mem_string_1} bytes")
print(f"Memory for the index variable (int i): {mem_index_variable} bytes")
print(f"Memory for the single character from the second string (int c): {mem_char_variable} bytes")
print(f"Total memory m = {mem_string_1} + {mem_index_variable} + {mem_char_variable} = {m} bytes")
