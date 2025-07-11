# The riddle describes a sequence of five numbers. I will determine the order
# based on the logical constraints provided in the riddle.

# The riddle states: "Number 1 is the first...". This places 1 at the start.
# "the fifth, who never had a thing and lastly follows." This places 5 at the end.
# "Number three is the strongest and it protects the last;". This places 3 just before the last number (5).
# So far, the sequence is: [1, ?, ?, 3, 5]

# The remaining numbers are 2 and 4. The clue is:
# "Number 4 likes two the most and it always goes behind."
# This means 4 comes immediately after 2. The only way to fit (2, 4) into the
# empty slots is to place 2 in the second position and 4 in the third.

# Therefore, the final ordered sequence is 1, 2, 4, 3, 5.

# The final code will print each number of this sequence.
final_sequence = [1, 2, 4, 3, 5]

# This part of the code fulfills the request to "output each number".
# It prints the members of the sequence separated by a space.
print(f"The correct order of the numbers is:")
print(*final_sequence)