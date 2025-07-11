# Plan:
# 1. Define the number of turns for each type of move in the sequence.
# 2. A pair of 180-degree turns (one away from the audience, one returning) is treated as one logical unit that shows the right side once, the left side once, and the back once.
# 3. A 360-degree turn (Barrel Turn) also shows each side and the back once.
# 4. Other moves in the sequence do not involve turning the body's orientation to the audience.
# 5. Sum the turns from all moves to get the final count.

# Number of turns from the Swivel Step Half Turn and Sunanda pair.
# This pair constitutes a turn away from and a turn back to the audience.
# This combination shows the right side once, the left side once, and the back once.
ssht_sunanda_pair_turns = {"right": 1, "left": 1, "back": 1}

# Number of turns from the Barrel Turn.
# This is a full 360-degree spin, showing each side and the back.
barrel_turn_turns = {"right": 1, "left": 1, "back": 1}

# Other moves (Swivel Step, Balancing Step, Figure 8) have no turns.
no_turns = {"right": 0, "left": 0, "back": 0}

# Calculate the total turns
total_right = ssht_sunanda_pair_turns["right"] + barrel_turn_turns["right"]
total_left = ssht_sunanda_pair_turns["left"] + barrel_turn_turns["left"]
total_back = ssht_sunanda_pair_turns["back"] + barrel_turn_turns["back"]

# Print the calculation for each orientation
print("Calculating the total turns:")
print(f"Right side turns = {ssht_sunanda_pair_turns['right']} (from Swivel/Sunanda pair) + {barrel_turn_turns['right']} (from Barrel Turn) = {total_right}")
print(f"Left side turns = {ssht_sunanda_pair_turns['left']} (from Swivel/Sunanda pair) + {barrel_turn_turns['left']} (from Barrel Turn) = {total_left}")
print(f"Back turns = {ssht_sunanda_pair_turns['back']} (from Swivel/Sunanda pair) + {barrel_turn_turns['back']} (from Barrel Turn) = {total_back}")

print(f"\nFinal Answer: {total_right}/{total_left}/{total_back}")