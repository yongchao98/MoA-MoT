def solve_poker_question():
  """
  This function identifies the worst suited jack to open from the button
  in a 100BB rake-free cash game based on GTO principles.
  """

  # In poker, the strength of a suited jack hand (JXs) is primarily
  # determined by its kicker (the 'X' card).
  # A higher kicker means higher equity and better straight possibilities.
  # The list of suited jacks from best to worst is:
  # JTs, J9s, J8s, J7s, J6s, J5s, J4s, J3s, J2s.

  # In a 100BB rake-free game, opening from the button is extremely profitable,
  # and a GTO (Game Theory Optimal) range includes all suited jacks.

  # The "worst" hand to open is the one with the lowest equity that is still
  # considered a standard, profitable raise.
  worst_suited_jack_to_open = "J2s"

  print(f"The worst suited jack you should open from the button is: {worst_suited_jack_to_open}")

solve_poker_question()