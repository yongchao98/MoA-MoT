def find_worst_suited_jack():
  """
  Determines the worst suited Jack to open from the button in a standard 100BB cash game.
  
  In poker, opening ranges from the button are very wide. All suited Jacks (Jxs)
  are generally considered profitable opens. The "worst" one is the hand with the
  lowest kicker, as it has the poorest playability and is most susceptible to
  domination.
  """
  
  hand_rank = "J"
  worst_kicker = "2"
  suit_modifier = "s" # 's' stands for suited
  
  # The abbreviated hand is formed by the rank, the kicker, and the suit modifier
  worst_hand = f"{hand_rank}{worst_kicker}{suit_modifier}"
  
  print(f"The worst suited Jack to open from the button is: {worst_hand}")
  
find_worst_suited_jack()