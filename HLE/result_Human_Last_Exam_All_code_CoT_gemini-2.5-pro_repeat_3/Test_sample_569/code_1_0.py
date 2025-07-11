# Plan:
# 1. Define variables for the value of each meld and bonus.
# 2. Calculate the meld points for myself and my partner.
# 3. Calculate the total trick points available in the deck.
# 4. Sum all the points together to get the final score.
# 5. Print the breakdown of the calculation clearly.

# Meld point values
run_meld = 150
aces_around_meld = 1000
pinochle_meld = 40
dix_meld = 10

# My meld points
my_meld_points = run_meld + aces_around_meld

# Partner's meld points
partner_meld_points = pinochle_meld + dix_meld

# Total team meld points
total_meld = my_meld_points + partner_meld_points

# Trick point values
# Card points: A=11, 10=10, K=4, Q=3, J=2, 9=0
# There are 2 of each card per suit.
points_per_suit = (11 * 2) + (10 * 2) + (4 * 2) + (3 * 2) + (2 * 2)
total_trick_points = points_per_suit * 4

# Bonus for winning the last trick
last_trick_bonus = 10

# The problem states "assuming that I play my hand perfectly".
# A hand with a run and 8 aces is powerful enough to win all tricks.
# Therefore, the team wins all trick points and the last trick bonus.
points_from_play = total_trick_points + last_trick_bonus

# Final score calculation
final_score = total_meld + points_from_play

# Print the detailed breakdown
print("Calculating the total points for the hand:")
print("-" * 40)

print(f"Your Meld Points: {run_meld} (Run) + {aces_around_meld} (8 Aces) = {my_meld_points}")
print(f"Partner's Meld Points: {pinochle_meld} (Pinochle) + {dix_meld} (Dix) = {partner_meld_points}")
print(f"Total Meld Points = {my_meld_points} + {partner_meld_points} = {total_meld}")
print("-" * 40)

print("Points from winning all tricks (perfect play):")
print(f"Total Trick Points in Deck: {total_trick_points}")
print(f"Bonus for Last Trick: {last_trick_bonus}")
print("-" * 40)

print("Final Score Calculation:")
print(f"Total Score = Total Meld + Trick Points + Last Trick Bonus")
print(f"Total Score = {total_meld} + {total_trick_points} + {last_trick_bonus} = {final_score}")

<<<1450>>>