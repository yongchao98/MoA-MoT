# Based on the riddle, I have deduced the correct order of the five natural numbers.
# The plan is to store these numbers in a list and then print them.

# The numbers involved are 0, 1, 2, 3, and 4.
# "Number three is the strongest" -> The value is 4.
# "Number 1 is the ... weakest" -> The value is 1.
# "the fifth, who never had a thing" -> The value is 0.
# "Number 2 is the best of them" -> The value is 3 (best of the remaining 2 and 3).
# "Number 4" -> The value is 2.

# Now, for the positions:
# "the fifth ... lastly follows" -> 0 is in the last position.
# "Number three ... protects the last" -> 4 is in the second to last position.
# "Number 1 is the first" -> 1 is in the first position.
# "Number 4 ... goes behind [Number 2]" -> The sequence 3, 2 fills the middle.

# This gives the final order.
final_sequence = [1, 3, 2, 4, 0]

# Printing the final sequence. The riddle asks to "Put the numbers in the correct order",
# and the instructions say to "output each number".
print("The correct order of the numbers is:")
# The format below prints each number separated by a comma and a space.
print(*final_sequence, sep=", ")