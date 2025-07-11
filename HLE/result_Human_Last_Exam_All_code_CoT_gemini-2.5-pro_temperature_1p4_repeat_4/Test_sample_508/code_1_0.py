import math

# Step 1 & 2: Define initial parameters and the observer's final position.
# Sources are defined as a list of tuples: (Name, Level @ 1m, x_pos, y_pos)
sources = [
    ("Dog", 55, -25, 0),
    ("Train", 110, 50, 0),
    ("Construction", 90, 0, 75),
    ("People", 75, 0, -10)
]
observer_final_pos = (0, 25)

print("Calculating the sound level at your new position (0, 25).\n")

# A dictionary to store the calculated sound level (in dB) for each source at the new location.
new_levels = {}

# Step 3 & 4: Calculate new distance and new sound level for each source.
print("First, let's find the sound level from each source at your new location:")
for name, level_at_1m, x_pos, y_pos in sources:
    # Calculate the distance from the observer's final position to the source.
    distance = math.sqrt((x_pos - observer_final_pos[0])**2 + (y_pos - observer_final_pos[1])**2)
    
    # Calculate the new sound level (L2) at the new distance (r2).
    # L2 = L1 - 20 * log10(r2 / r1). Since r1 = 1m, this simplifies to L1 - 20 * log10(r2).
    level_at_new_distance = level_at_1m - 20 * math.log10(distance)
    
    # Store the result for the final calculation.
    new_levels[name] = level_at_new_distance
    
    print(f"- The {name} is {distance:.2f} meters away, creating a sound level of {level_at_new_distance:.2f} dB.")

# Step 5 & 6: Formulate the equation to sum the intensities and calculate the final result.
# We will use the calculated levels to build the final equation string.
# The formula is L_total = 10 * log10(10^(L_dog/10) + 10^(L_train/10) + ...)
sum_of_intensities_terms = []
for name in new_levels:
    term = f"10^({new_levels[name]:.2f}/10)"
    sum_of_intensities_terms.append(term)

equation_str = f"Total dB = 10 * log10({ ' + '.join(sum_of_intensities_terms) })"

# Calculate the actual total intensity and then the total dB level.
total_intensity = sum(10**(level / 10) for level in new_levels.values())
total_db = 10 * math.log10(total_intensity)

print("\nTo find the total sound level, we convert each dB value to intensity, sum them, and convert back to dB.")
print("The final equation is:")
print(equation_str)

print(f"\nThe total sound level you hear is: {total_db:.2f} dB")

# Final answer in the specified format
print(f"\n<<<{total_db:.2f}>>>")