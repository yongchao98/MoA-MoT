# Step 1: Define the classification for the first nanoribbon.
# Edge is Armchair ('A'), Width (N) is 8, Band type is semiconducting ('1').
edge_1 = 'A'
width_1 = 8
band_1 = 1

# Step 2: Define the classification for the second nanoribbon.
# Edge is Zigzag ('Z'), Width (N) is 5, Band type is metallic ('0').
edge_2 = 'Z'
width_2 = 5
band_2 = 0

# Step 3: Define the classification for the third nanoribbon.
# Edge is Zigzag ('Z'), Width (N) is 7, Band type is metallic ('0').
edge_3 = 'Z'
width_3 = 7
band_3 = 0

# Step 4: Construct the individual classification strings.
class_1 = f"{edge_1}{width_1}{band_1}"
class_2 = f"{edge_2}{width_2}{band_2}"
class_3 = f"{edge_3}{width_3}{band_3}"

# Step 5: Concatenate the strings and print the final result.
final_classification = class_1 + class_2 + class_3
print(final_classification)