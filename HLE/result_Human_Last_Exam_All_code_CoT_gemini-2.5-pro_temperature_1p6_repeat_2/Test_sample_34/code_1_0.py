# The hand in question is a suited Jack.
# We are looking for the weakest one that is still a standard open-raise from the button.
# According to GTO (Game Theory Optimal) charts for a 100BB rake-free cash game,
# the button's opening range is very wide and includes all suited Jacks.

# The highest card is a Jack.
high_card = "J"

# The lowest possible kicker that makes the hand a profitable open-raise is a 2.
low_card = "2"

# The "s" indicates that the two cards are of the same suit.
suit_info = "s"

# Combining these components gives us the abbreviated hand name.
worst_suited_jack_to_open = high_card + low_card + suit_info

print("The worst suited jack you should open from the button is:")
print(worst_suited_jack_to_open)