import math

# Step 1: Explain the interpretation of the clue "Кома" (Koma) as "comma",
# suggesting coordinates derived from mathematical constants pi and e.
print("The clue 'Кома' (Koma) is interpreted as a pun on 'comma', indicating a location defined by coordinates.")
print("These coordinates can be calculated using the fundamental constants π and e, which a software engineer would know.\n")

# Step 2: Define the equations for Latitude and Longitude.
print("--- The Equation for the Location ---")

# Step 3: Calculate and print the Latitude value.
# The formula is based on the conversion of 1 radian to degrees.
lat_numerator = 180
pi_val = math.pi
latitude = lat_numerator / pi_val

print(f"Latitude = {lat_numerator} / π")
print(f"Latitude = {lat_numerator} / {pi_val}")
print(f"Result: {latitude:.4f}° N\n")

# Step 4: Calculate and print the Longitude value.
# The formula uses Euler's number (e) raised to a simple integer power.
e_val = math.e
lon_power = 4
longitude = e_val ** lon_power

print(f"Longitude = e ^ {lon_power}")
print(f"Longitude = {e_val} ^ {lon_power}")
print(f"Result: {longitude:.4f}° E\n")

# Step 5: State the final conclusion.
print(f"The calculated coordinates are approximately {latitude:.2f}° N, {longitude:.2f}° E.")
print("This point is located in Perm Krai.")
