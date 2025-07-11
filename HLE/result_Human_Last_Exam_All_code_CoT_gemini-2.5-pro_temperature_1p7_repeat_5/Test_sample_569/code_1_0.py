# Step 1: Define and calculate your meld points.
my_run_points = 15
my_eight_aces_points = 100
my_total_meld = my_run_points + my_eight_aces_points

# Step 2: Define and calculate your partner's meld points.
# A single pinochle is worth 4 points.
partner_pinochle_points = 4
# Each 9 of trump (dix) is worth 1 point. Your partner has two.
partner_dix_points = 2
partner_total_meld = partner_pinochle_points + partner_dix_points

# Step 3: Calculate the team's total meld points.
team_meld_total = my_total_meld + partner_total_meld

# Step 4: Calculate the trick points.
# A 48-card Pinochle deck has two of each card from 9 to Ace for each suit.
# This means there are 8 Aces, 8 Tens, and 8 Kings in total.
# Your hand (run + 8 aces) gives you an unbeatable trump suit and the top card in every other suit.
# With perfect play, you will win all 12 tricks. This means your team captures all point-bearing cards.
ace_trick_points = 8
ten_trick_points = 8
king_trick_points = 8
last_trick_bonus = 1
team_trick_total = ace_trick_points + ten_trick_points + king_trick_points + last_trick_bonus

# Step 5: Calculate the grand total for the hand.
grand_total = team_meld_total + team_trick_total

# Final Output
# We will print the breakdown of the final score.
print("Here is the breakdown of the points for the hand:")
print("\nMeld Points:")
print(f"Your Run: {my_run_points}")
print(f"Your 8 Aces: {my_eight_aces_points}")
print(f"Partner's Pinochle: {partner_pinochle_points}")
print(f"Partner's 2 Nines of Trump: {partner_dix_points}")

print("\nTrick Points (from winning all tricks):")
print(f"Points from Aces: {ace_trick_points}")
print(f"Points from 10s: {ten_trick_points}")
print(f"Points from Kings: {king_trick_points}")
print(f"Bonus for Last Trick: {last_trick_bonus}")

print("\nFinal Calculation:")
# Print the full equation as requested
print(f"{my_run_points} (Run) + {my_eight_aces_points} (8 Aces) + {partner_pinochle_points} (Pinochle) + {partner_dix_points} (2 Dixes) + {ace_trick_points} (Aces) + {ten_trick_points} (10s) + {king_trick_points} (Kings) + {last_trick_bonus} (Last Trick) = {grand_total}")
print(f"\nYour team will earn a total of {grand_total} points this hand.")
