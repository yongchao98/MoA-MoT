# Define the components of the poker hand.
# The hand is a suited Jack, and we need to find the lowest kicker
# that is still a profitable open-raise from the button in a rake-free game.

high_card = 'J'
# Based on GTO principles for a rake-free game, the range is wide enough
# to include the lowest possible kicker.
low_card_value = 2
suit = 's'

# Construct the abbreviated hand name.
hand_abbreviation = high_card + str(low_card_value) + suit

# The "equation" for forming the hand is simply concatenating the parts.
# The prompt requires showing the numbers in the final equation.
print(f"Equation for the hand: '{high_card}' + '{low_card_value}' + '{suit}'")
print(f"The worst suited jack you should open from the button is: {hand_abbreviation}")