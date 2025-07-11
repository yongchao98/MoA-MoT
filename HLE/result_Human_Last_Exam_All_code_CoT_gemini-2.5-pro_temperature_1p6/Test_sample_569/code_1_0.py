# Plan:
# 1. Calculate meld points for both partners.
# 2. Calculate the total points available in the deck from trick-taking.
# 3. Add the last trick bonus.
# 4. Sum all points for the final score.

# Step 1: Calculate Meld Points

# Your meld points
my_run = 15  # A family (or run) in trump is 15 points.
my_aces = 100 # All 8 aces (100 aces) is 100 points.

# Partner's meld points
partner_pinochle = 4 # A pinochle (Jack of Diamonds, Queen of Spades) is 4 points.
partner_dix = 10 # A nine of trump (dix) is worth 10 points.

total_meld_points = my_run + my_aces + partner_pinochle + partner_dix

# Step 2: Calculate Total Trick Points

# Card values for tricks
ace_value = 11
ten_value = 10
king_value = 4
queen_value = 3
jack_value = 2

# There are 8 of each card in a 48-card Pinochle deck.
num_of_each_card = 8

total_card_points = num_of_each_card * (ace_value + ten_value + king_value + queen_value + jack_value)

# Step 3: Add Last Trick Bonus
# With your hand, perfect play means winning all tricks, including the last one.
last_trick_bonus = 10

total_trick_points = total_card_points + last_trick_bonus

# Step 4: Calculate Final Score
final_score = total_meld_points + total_trick_points

# The problem states "other players play a 9 and two jacks of trump."
# In a 48-card deck, there are only two jacks of trump. If you hold one for your run,
# it is impossible for the other players to play two.
# However, this detail does not change the outcome, as "perfect play" with your
# described hand means you will win all 12 tricks regardless of how the opponents play.
# We proceed by calculating the score based on winning all possible points.

print("To calculate your total score, we sum the points from melds and tricks.")
print("\nYour Team's Meld Points:")
print(f"  Your run: {my_run}")
print(f"  Your 100 aces: {my_aces}")
print(f"  Partner's pinochle: {partner_pinochle}")
print(f"  Partner's dix (9 of trump): {partner_dix}")

print("\nPoints from Winning All Tricks:")
print(f"  Total points from cards in the deck: {total_card_points}")
print(f"  Bonus for winning the last trick: {last_trick_bonus}")

print("\nFinal Score Calculation:")
print(f"{my_run} (run) + {my_aces} (aces) + {partner_pinochle} (pinochle) + {partner_dix} (dix) + {total_card_points} (cards) + {last_trick_bonus} (last trick) = {final_score}")
<<<379>>>