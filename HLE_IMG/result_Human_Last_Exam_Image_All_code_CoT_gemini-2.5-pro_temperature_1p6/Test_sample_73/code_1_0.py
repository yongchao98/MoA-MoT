def get_stereochemical_assignments():
  """
  This function provides the stereochemical assignments (R/S) for the four
  stereocenters in the reaction scheme, moving from left to right.
  """
  # Assignment for the first stereocenter (acyl chloride)
  stereocenter1 = "R"

  # Assignment for the second stereocenter (alcohol)
  stereocenter2 = "S"

  # Assignment for the third stereocenter (product, from alcohol)
  stereocenter3 = "S"

  # Assignment for the fourth stereocenter (product, from acyl chloride)
  stereocenter4 = "S"

  # Print the final sequence
  print(f"{stereocenter1}, {stereocenter2}, {stereocenter3}, {stereocenter4}")

get_stereochemical_assignments()
<<<R, S, S, S>>>