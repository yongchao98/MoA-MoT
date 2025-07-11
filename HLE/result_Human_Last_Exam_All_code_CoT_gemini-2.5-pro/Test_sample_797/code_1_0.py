# This script calculates the number of turns based on a sequence of FCBD® dance moves.

# Step 1: Define the turns for each move in the sequence.
# (Right turns, Left turns, Back turns)
# Note: Some moves are interpreted as showing side profiles without showing the back.
swivel_step_half_turn = (1, 1, 0)
sunanda = (1, 1, 1)
balancing_step = (1, 1, 0)
barrel_turn = (1, 1, 1)
# Swivel Step and Figure 8 have no turns, so they contribute 0 and are omitted from the equation.

# Step 2: Calculate the total turns for each direction.
total_right_turns = swivel_step_half_turn[0] + sunanda[0] + balancing_step[0] + barrel_turn[0]
total_left_turns = swivel_step_half_turn[1] + sunanda[1] + balancing_step[1] + barrel_turn[1]
total_back_turns = swivel_step_half_turn[2] + sunanda[2] + balancing_step[2] + barrel_turn[2]

# Step 3: Print the breakdown of the calculation and the final result.
print("Calculating total turns for Right, Left, and Back.")

# Print the equation for Right turns
print("\nTotal Right Turns:")
print(f"  {swivel_step_half_turn[0]} (from Swivel Step Half Turn)")
print(f"+ {sunanda[0]} (from Sunanda)")
print(f"+ {balancing_step[0]} (from Balancing Step)")
print(f"+ {barrel_turn[0]} (from Barrel Turn)")
print(f"--------------------")
print(f"= {total_right_turns}")

# Print the equation for Left turns
print("\nTotal Left Turns:")
print(f"  {swivel_step_half_turn[1]} (from Swivel Step Half Turn)")
print(f"+ {sunanda[1]} (from Sunanda)")
print(f"+ {balancing_step[1]} (from Balancing Step)")
print(f"+ {barrel_turn[1]} (from Barrel Turn)")
print(f"--------------------")
print(f"= {total_left_turns}")

# Print the equation for Back turns
print("\nTotal Back Turns:")
print(f"  {swivel_step_half_turn[2]} (from Swivel Step Half Turn)")
print(f"+ {sunanda[2]} (from Sunanda)")
print(f"+ {balancing_step[2]} (from Balancing Step)")
print(f"+ {barrel_turn[2]} (from Barrel Turn)")
print(f"--------------------")
print(f"= {total_back_turns}")


print(f"\nFinal Answer: The dancer turns her/his right side {total_right_turns} times, left side {total_left_turns} times, and back {total_back_turns} times.")