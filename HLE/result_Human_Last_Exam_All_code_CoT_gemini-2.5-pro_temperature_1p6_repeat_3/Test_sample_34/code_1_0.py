def find_worst_suited_jack_open():
  """
  This function determines the worst suited Jack to open from the button
  in a 100BB rake-free cash game based on standard GTO strategy.

  In poker GTO charts for the button position, the opening range is very wide.
  For suited Jacks (JXs), all combinations from JQs down to J2s are
  typically considered profitable opens.

  The strength of these hands decreases as the kicker (the 'X') gets lower.
  Therefore, the hand with the lowest EV (Expected Value) that is still a
  standard open-raise is the one with the lowest kicker.
  """
  worst_hand = "J2s"
  print(worst_hand)

find_worst_suited_jack_open()