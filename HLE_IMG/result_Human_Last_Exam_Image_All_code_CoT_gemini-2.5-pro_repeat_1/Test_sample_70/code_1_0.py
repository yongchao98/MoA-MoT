def identify_pericyclic_reactions():
  """
  This function explains the two pericyclic reactions involved in the given
  thermal transformation of bicyclo[6.2.0]decatetraene to cis-9,10-dihydronaphthalene.
  """
  
  # Explanation of the first reaction
  reaction_1_name = "Electrocyclic Ring Opening"
  reaction_1_description = ("The four-membered cyclobutene ring of the starting material "
                            "undergoes a thermal 4π-electron electrocyclic ring opening. "
                            "This reaction forms a highly unstable 10-membered ring intermediate "
                            "(cyclodeca-1,3,5,7,9-pentaene).")

  # Explanation of the second reaction
  reaction_2_name = "Electrocyclic Ring Closing"
  reaction_2_description = ("The intermediate, which contains a 1,3,5-hexatriene system, "
                            "immediately undergoes a thermal 6π-electron electrocyclic ring closing. "
                            "This second reaction forms the final, more stable product, "
                            "cis-9,10-dihydronaphthalene.")

  # Print the results
  print("The thermal transformation proceeds via a sequence of two specific pericyclic reactions:\n")
  print(f"1. {reaction_1_name}:")
  print(f"   {reaction_1_description}\n")
  print(f"2. {reaction_2_name}:")
  print(f"   {reaction_2_description}")

identify_pericyclic_reactions()