def get_sound_hierarchy():
  """
  Determines and prints the hierarchy of sound importance in a dance scene.
  
  In a dance scene, music is the foundational element that drives the action and emotion.
  Sound effects provide texture and realism (e.g., footsteps, rustling clothes).
  Speech is typically the least important element as the focus is on movement and music.
  """
  
  # The hierarchy is Music > Sound Effects > Speech.
  # We list them in this order.
  hierarchy_list = ["Music", "Sound Effects", "Speech"]
  
  # Format the list into a string separated by a comma and a space.
  result_string = ", ".join(hierarchy_list)
  
  # Print the final formatted string.
  print(result_string)

get_sound_hierarchy()