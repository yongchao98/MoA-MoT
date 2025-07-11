def get_hook_flash_champions():
  """
  This function identifies and prints the names of champions from the given list
  who can perform a 'hook-flash' in League of Legends.

  A hook-flash is a mechanic where a champion uses Flash during the cast
  animation of their hook ability to extend its range or change its angle.
  """
  
  # The list of champions to check from.
  champion_pool = ["Blitzcrank", "Nautilus", "Pyke", "Thresh"]

  # All champions in the provided list are capable of performing a hook-flash.
  # Blitzcrank (Q - Rocket Grab)
  # Nautilus (Q - Dredge Line)
  # Pyke (Q - Bone Skewer)
  # Thresh (Q - Death Sentence)
  capable_champions = ["Blitzcrank", "Nautilus", "Pyke", "Thresh"]

  # Join the names with a comma for the final output.
  result = ",".join(capable_champions)
  
  print(result)

get_hook_flash_champions()