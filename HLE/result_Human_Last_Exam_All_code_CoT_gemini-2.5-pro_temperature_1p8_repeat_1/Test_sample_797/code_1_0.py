# Define the number of turns (Right, Left, Back) for each FCBD® move.
# A full 360-degree turn presents the right side, back, and left side once.
# A torso twist is assumed to present the right and left sides.
moves_definitions = {
    "Swivel Step": (0, 0, 0),
    "Swivel Step Half Turn": (1, 1, 1),
    "Sunanda": (1, 1, 1),
    "Balancing Step": (1, 1, 0),
    "Figure 8": (0, 0, 0),
    "Barrel Turn": (1, 1, 1)
}

# Define the sequence from the problem description.
# Note: Repetitions for Figure 8 are handled in the calculation.
sequence_moves = [
    "Swivel Step",
    "Swivel Step Half Turn",
    "Sunanda",
    "Balancing Step",
    "Figure 8",
    "Barrel Turn"
]

# Get the turn counts for each move in the sequence
ss_r, ss_l, ss_b = moves_definitions["Swivel Step"]
ssht_r, ssht_l, ssht_b = moves_definitions["Swivel Step Half Turn"]
sunanda_r, sunanda_l, sunanda_b = moves_definitions["Sunanda"]
balancing_r, balancing_l, balancing_b = moves_definitions["Balancing Step"]
figure8_r, figure8_l, figure8_b = moves_definitions["Figure 8"] # For 8 repetitions, this is still 0
barrel_r, barrel_l, barrel_b = moves_definitions["Barrel Turn"]

# Calculate the total turns for each orientation
total_right_turns = ss_r + ssht_r + sunanda_r + balancing_r + (figure8_r * 8) + barrel_r
total_left_turns = ss_l + ssht_l + sunanda_l + balancing_l + (figure8_l * 8) + barrel_l
total_back_turns = ss_b + ssht_b + sunanda_b + balancing_b + (figure8_b * 8) + barrel_b

# Print the final calculations and results, showing each number in the equation.
print("Calculating the total number of turns...\n")

print("Total Right Side Turns:")
print(f"{ss_r} + {ssht_r} + {sunanda_r} + {balancing_r} + {figure8_r * 8} + {barrel_r} = {total_right_turns}")

print("\nTotal Left Side Turns:")
print(f"{ss_l} + {ssht_l} + {sunanda_l} + {balancing_l} + {figure8_l * 8} + {barrel_l} = {total_left_turns}")

print("\nTotal Back to Audience Turns:")
print(f"{ss_b} + {ssht_b} + {sunanda_b} + {balancing_b} + {figure8_b * 8} + {barrel_b} = {total_back_turns}")

print(f"\nFinal Answer: {total_right_turns}/{total_left_turns}/{total_back_turns}")