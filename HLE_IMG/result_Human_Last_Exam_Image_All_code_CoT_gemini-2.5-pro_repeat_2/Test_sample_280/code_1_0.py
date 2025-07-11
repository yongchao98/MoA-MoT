# The provided image is a well-known Piet program. For it to work correctly,
# the interpreter must be set to use a "codel" size of 16 pixels.
# The program's logic is to arithmetically construct the ASCII character codes
# for the string "hello" and then print them to standard output.

h_ascii = 104
e_ascii = 101
l_ascii = 108
o_ascii = 111

# The "equation" to form the word "hello" is based on these ASCII values.
# As requested, we will output each number used.
print(f"{h_ascii}")
print(f"{e_ascii}")
print(f"{l_ascii}")
print(f"{l_ascii}")
print(f"{o_ascii}")

# The final, concatenated output of the Piet program is the string itself.
final_output = chr(h_ascii) + chr(e_ascii) + chr(l_ascii) + chr(l_ascii) + chr(o_ascii)
print(final_output)
