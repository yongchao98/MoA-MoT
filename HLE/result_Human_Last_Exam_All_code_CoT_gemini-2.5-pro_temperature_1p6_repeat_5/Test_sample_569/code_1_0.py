# Plan: Calculate total score from a Pinochle hand.
# The total score is the sum of meld points and trick points.

# 1. Calculate meld points for our team.
# My meld:
# - A run (A, K, Q, J, 10 of trump) is worth 15 points.
my_run_meld = 15
# - Having all 8 aces (two of each suit) is "Double Aces Around" and is worth 100 points.
my_aces_meld = 100
my_total_meld = my_run_meld + my_aces_meld

# Partner's meld:
# - A pinochle (Queen of Spades and Jack of Diamonds) is worth 4 points.
partner_pinochle_meld = 4
# - Each 9 of trump (called a 'dix') is worth 1 point. Partner has two.
partner_dix_meld = 2 * 1
partner_total_meld = partner_pinochle_meld + partner_dix_meld

# Total team meld points
total_meld = my_total_meld + partner_total_meld

# 2. Calculate trick points.
# With all 8 aces and trump superiority, perfect play guarantees winning all 12 tricks.
# This means we capture all "counter" cards.
# There are 8 Aces, 8 Tens, and 8 Kings in the deck, each worth 1 point.
num_aces = 8
num_tens = 8
num_kings = 8
trick_counters = num_aces + num_tens + num_kings

# There is a 1-point bonus for winning the last trick.
last_trick_bonus = 1

# Total trick points
total_trick_points = trick_counters + last_trick_bonus

# 3. Calculate the total score for the hand.
total_score = total_meld + total_trick_points

# Print the final breakdown of the score.
# The user wants each number in the final equation.
print(f"Your team's score is calculated by adding all meld and trick points:")
print(f"{my_run_meld} (My Run) + {my_aces_meld} (My 8 Aces) + {partner_pinochle_meld} (Partner's Pinochle) + {partner_dix_meld} (Partner's 2 Nines) + {trick_counters} (All Counters) + {last_trick_bonus} (Last Trick) = {total_score}")
