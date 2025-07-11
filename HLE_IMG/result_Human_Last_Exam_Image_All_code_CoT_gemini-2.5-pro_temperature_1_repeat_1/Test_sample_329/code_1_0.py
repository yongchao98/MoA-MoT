# Based on the analysis of the four simulations, here are the identified unique parameters and their codes.

# Simulation 1: The unique parameter is u1 (wave number of the left slit).
sim1_code = 4

# Simulation 2: The unique parameter is u3 (wave number of the right slit).
sim2_code = 6

# Simulation 3: The unique parameter is h2 (height/amplitude of the middle slit).
sim3_code = 8

# Simulation 4: The unique parameter is h2 (height/amplitude of the middle slit).
sim4_code = 8

# The problem asks for the four-digit sequence.
print(f"The code for Simulation 1 is: {sim1_code}")
print(f"The code for Simulation 2 is: {sim2_code}")
print(f"The code for Simulation 3 is: {sim3_code}")
print(f"The code for Simulation 4 is: {sim4_code}")

final_sequence = f"{sim1_code}{sim2_code}{sim3_code}{sim4_code}"
print(f"The final four-digit sequence is: {final_sequence}")
