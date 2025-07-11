# My meld components
my_run_points = 15
my_aces_points = 100

# Partner's meld components
partner_pinochle_points = 4
# Two 9s of trump (dix), each worth 10 points
partner_dix_points = 2 * 10

# Trick point components
# In a 48-card deck, there are 8 aces, 8 tens, and 8 kings.
# Assuming a perfect hand, the partnership wins all of them.
# Each card is worth 1 point.
trick_card_points = 8 + 8 + 8
last_trick_bonus = 1

# Calculate the total score by summing all components
total_score = (
    my_run_points
    + my_aces_points
    + partner_pinochle_points
    + partner_dix_points
    + trick_card_points
    + last_trick_bonus
)

# Print the final equation with all its components
print("The total points are calculated by summing the points from your meld, your partner's meld, and the tricks.")
print("Final Equation:")
print(
    f"{my_run_points} (run) + {my_aces_points} (8 aces) + "
    f"{partner_pinochle_points} (pinochle) + {partner_dix_points} (2 dix) + "
    f"{trick_card_points} (card points) + {last_trick_bonus} (last trick) = {total_score}"
)
print(f"Total points earned this hand: {total_score}")
