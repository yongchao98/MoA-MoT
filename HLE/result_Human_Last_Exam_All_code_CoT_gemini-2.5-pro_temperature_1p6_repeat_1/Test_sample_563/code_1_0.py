def get_automorphism_group_counts():
  """
  Provides the number of isomorphism classes of automorphism groups for
  compact, connected Riemann surfaces of genus g=2, 3, and 4.
  """
  # Number of groups for genus g=2 is 12.
  genus_2_groups = 12

  # Number of groups for genus g=3 is 15.
  genus_3_groups = 15

  # Number of groups for genus g=4 is 23.
  genus_4_groups = 23

  result = [genus_2_groups, genus_3_groups, genus_4_groups]
  print(result)

get_automorphism_group_counts()