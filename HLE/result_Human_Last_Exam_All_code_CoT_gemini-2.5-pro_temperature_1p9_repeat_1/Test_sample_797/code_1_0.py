import collections

# Plan:
# We will analyze each FCBD® move for the number of right-side, left-side,
# and back turns it contributes. A "turn" is counted as an action of presenting
# a side or the back to the audience, based on the definition of the move.

# Initialize turn counters
total_turns = collections.Counter()

# --- Breakdown of the FCBD® Sequence ---

print("Analyzing the dance sequence step-by-step:\n")

# Move 1: Swivel Step (8 counts)
# This move has no turns.
move1_turns = {'right': 0, 'left': 0, 'back': 0}

# Move 2: Swivel Step Half Turn (16 counts)
# A 180-degree turn ending with the back to the audience (+1 back).
# The turn itself involves passing one side. As the direction is not specified,
# we count it abstractly as 1 right and 1 left turn.
move2_turns = {'right': 1, 'left': 1, 'back': 1}
total_turns.update(move2_turns)
print(f"Swivel Step Half Turn adds: {move2_turns['right']} right, {move2_turns['left']} left, {move2_turns['back']} back")

# Move 3: Sunanda (once)
# Consists of two 180-degree spins. The spins over each shoulder present both
# profiles (+1 right, +1 left). The second spin lands facing back (+1 back).
move3_turns = {'right': 1, 'left': 1, 'back': 1}
total_turns.update(move3_turns)
print(f"Sunanda adds:               {move3_turns['right']} right, {move3_turns['left']} left, {move3_turns['back']} back")

# Move 4: Balancing Step (once)
# This move has no turns.
move4_turns = {'right': 0, 'left': 0, 'back': 0}

# Move 5: Figure 8 (8 times)
# This move has no turns.
move5_turns = {'right': 0, 'left': 0, 'back': 0}

# Move 6: Barrel Turn (once)
# A 360-degree spin that presents both side profiles (+1 right, +1 left).
# It ends in the starting orientation, so no new back turn is counted.
move6_turns = {'right': 1, 'left': 1, 'back': 0}
total_turns.update(move6_turns)
print(f"Barrel Turn adds:           {move6_turns['right']} right, {move6_turns['left']} left, {move6_turns['back']} back")

# --- Final Calculation ---
print("\n---")
print("Summing the contributions for the final count:\n")

# Print the final equation for each turn type
print(f"Total Right Turns = {move2_turns['right']} (from Swivel Step Half Turn) + {move3_turns['right']} (from Sunanda) + {move6_turns['right']} (from Barrel Turn) = {total_turns['right']}")
print(f"Total Left Turns  = {move2_turns['left']} (from Swivel Step Half Turn) + {move3_turns['left']} (from Sunanda) + {move6_turns['left']} (from Barrel Turn) = {total_turns['left']}")
print(f"Total Back Turns  = {move2_turns['back']} (from Swivel Step Half Turn) + {move3_turns['back']} (from Sunanda) + {move6_turns['back']} (from Barrel Turn) = {total_turns['back']}")

# --- Final Answer ---
print("\n---")
print(f"Final Answer: The dancer turns her/his right side/left side/back to the audience {total_turns['right']}/{total_turns['left']}/{total_turns['back']} times.")
