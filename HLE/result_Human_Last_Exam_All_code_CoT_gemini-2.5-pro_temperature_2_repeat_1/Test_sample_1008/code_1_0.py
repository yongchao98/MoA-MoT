def get_sound_hierarchy():
  """
  This function determines and prints the hierarchy of sound importance
  for a typical dance scene in the film industry.

  The ranking is based on the foundational role each sound element plays in constructing the scene.
  - Music: The primary driver, setting the rhythm and mood.
  - Sound Effects: Adds realism, texture, and physical presence (footsteps, clothing sounds).
  - Speech: Provides narrative but is often secondary to the music and action.
  - Silence: Used for dramatic effect, but is the absence of sound and thus the least foundational element.
  """
  # The list is ordered from most to least important.
  sound_elements_hierarchy = ["Music", "Sound Effects", "Speech", "Silence"]

  # Format the list into the required comma-separated string.
  output_string = ", ".join(sound_elements_hierarchy)

  # Print the final ranked order.
  print(output_string)

get_sound_hierarchy()