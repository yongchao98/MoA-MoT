def azobenzene_table_action(time_of_day):
  """
  Describes the action of an azobenzene-shaped table at different times of day.
  """
  if time_of_day.lower() == "sunrise":
    print("Action for Sunrise:")
    print("The sun's UV light hits the table.")
    print("The molecule isomerizes from the stable 'trans' form to the 'cis' form.")
    print("Result: The table bends from a straight shape into a V-shape.")
    print("\nVisual Representation:")
    print("trans (straight): \\_N=N_/   ----UV Light---->   cis (bent): \\_N")
    print("                                                     //")
    print("                                                    N__/")

  elif time_of_day.lower() == "sunset":
    print("Action for Sunset:")
    print("The sun's UV light disappears.")
    print("The molecule thermally relaxes from the 'cis' form back to the stable 'trans' form.")
    print("Result: The table straightens out from a V-shape back to its elongated shape.")
    print("\nVisual Representation:")
    print("cis (bent): \\_N   ----Darkness/Heat---->   trans (straight): \\_N=N_/")
    print("           //")
    print("          N__/")
  else:
    print("Please specify 'sunrise' or 'sunset'.")

# Demonstrate the actions
azobenzene_table_action("sunrise")
print("\n" + "="*40 + "\n")
azobenzene_table_action("sunset")