#
# Pinochle Hand Score Calculation
#

# --- Meld Point Calculation ---

# My meld points
my_run_points = 15  # A, 10, K, Q, J of trump
my_aces_points = 100 # All 8 Aces (double aces around)
my_total_meld = my_run_points + my_aces_points

# Partner's meld points
partner_pinochle_points = 4 # Q of Spades and J of Diamonds
partner_nines_points = 20 # Both 9s of trump (2 * 10 points)
partner_total_meld = partner_pinochle_points + partner_nines_points

# --- Trick Point Calculation ---
# Assuming perfect play with a dominant hand, the team will win all tricks.
# There are 8 Aces, 8 Tens, and 8 Kings in the deck, each worth 1 point.
counter_card_points = 24
last_trick_bonus = 1
total_trick_points = counter_card_points + last_trick_bonus

# --- Final Score Calculation ---
# The total score is the sum of all meld and all trick points for the team.
final_score = my_total_meld + partner_total_meld + total_trick_points

print("Your total meld points: {}".format(my_total_meld))
print("Your partner's total meld points: {}".format(partner_total_meld))
print("Total points from tricks: {}".format(total_trick_points))
print("\nThe final score is calculated as follows:")
print("{} (Your Meld) + {} (Partner's Meld) + {} (Tricks) = {}".format(my_total_meld, partner_total_meld, total_trick_points, final_score))

print("\nFinal Answer:")
print("The team will earn {} points this hand.".format(final_score))