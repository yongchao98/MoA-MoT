# Step 1: Define the numbers based on the riddle's logic.
# "Number 1 is the first to have something..."
# The sequence is 5 numbers long. Let's call them p1, p2, p3, p4, p5.
p1 = 1

# "...the fifth, who never had a thing and lastly follows."
p5 = 5

# "Number three is the strongest and it protects the last;"
# This means number 3 is in the 4th position, just before the last number.
p4 = 3

# "Number 4 likes two the most and it always goes behind."
# This means the pair 2, 4 must appear in the sequence.
# The only available spots are the 2nd and 3rd positions.
p2 = 2
p3 = 4

# Step 2: Print the final ordered sequence.
# The riddle asks to "Put the numbers in the correct order".
print("The correct order of the numbers is:")
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)