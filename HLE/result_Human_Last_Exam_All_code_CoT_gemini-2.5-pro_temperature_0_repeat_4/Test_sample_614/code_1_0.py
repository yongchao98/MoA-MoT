# The riddle asks to order five natural numbers based on a set of clues.
# Let's determine the numbers and their positions.

# Clue 5: "the fifth, who never had a thing and lastly follows."
# "Never had a thing" logically points to the number 0.
# "Lastly follows" places it at the end of the sequence.
fifth_number = 0

# Clue 3: "Number three ... protects the last;"
# This means the number 3 is in the fourth position, right before the last number (0).
fourth_number = 3

# Clue 1: "Number 1 is the first to have something..."
# This places the number 1 at the beginning of the sequence.
first_number = 1

# At this point, the sequence looks like: 1, _, _, 3, 0.
# The numbers left to place are 2 and 4. The open slots are the second and third positions.

# Clue 4: "Number 4 likes two the most and it always goes behind."
# This establishes the order for 2 and 4: 2 must come before 4.
# Therefore, the second number is 2 and the third number is 4.
second_number = 2
third_number = 4

# The complete sequence is 1, 2, 4, 3, 0.

# Let's check the remaining parts of the clues for consistency.
# - "Number 1 ... is the weakest": 1 is numerically the weakest (excluding 0). This fits.
# - "Number three is the strongest": This seems contradictory since 4 > 3. However, the "bloodline" clue helps.
#   The "bloodline" can be interpreted as Roman numerals (I, II, III, IV), which rule the numbers 1, 2, 3, and 4.
#   In this context, 3 (III) is the "strongest" number made by simple addition of 'I's before needing a new symbol or subtraction. This resolves the conflict.
# - "Number 2 is the best of them": This is subjective, but 2 is the only even prime, a unique quality that can be considered "best".

# The final order is derived from the positional clues.
print("The final sequence is:")
print(first_number, second_number, third_number, fourth_number, fifth_number)