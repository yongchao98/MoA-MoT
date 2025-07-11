def identify_histopathology(image_features):
  """
  This function provides a diagnosis based on key features observed in a histopathological image.

  Args:
    image_features: A list of strings describing the observed features.

  Returns:
    A string containing the diagnosis.
  """
  if "cellular crescent in glomerulus" in image_features:
    return "Crescentic Glomerulonephritis"
  else:
    return "Diagnosis cannot be determined from the provided features."

# The provided image prominently displays a hypercellular crescent filling Bowman's space
# and compressing the glomerular tuft. This is the pathognomonic feature.
features_observed = ["cellular crescent in glomerulus", "glomerular tuft compression", "tubular injury"]
diagnosis = identify_histopathology(features_observed)

print("Histopathological Diagnosis:")
print(diagnosis)