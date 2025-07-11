# Plan:
# 1. Define variables for each component of the meld score for both players.
# 2. Define variables for the points available in the tricks (counters and last trick bonus).
# 3. Sum all these components to get the total score.
# 4. Print the detailed equation and the final result.

# --- Meld Point Components ---
# Your meld
my_run = 15
my_aces = 100

# Partner's meld
partner_pinochle = 4
partner_nines = 2

# --- Trick Point Components ---
# The phrase "play my hand perfectly" with this strong hand implies taking all tricks.
# There are 8 Aces, 8 Tens, and 8 Kings in the deck, each worth 1 point.
all_counters = 24
last_trick_bonus = 1

# --- Total Score Calculation ---
total_score = my_run + my_aces + partner_pinochle + partner_nines + all_counters + last_trick_bonus

# --- Final Output ---
print("The total number of points earned is calculated as follows:")
print(f"{my_run} (your run) + {my_aces} (your 8 aces) + {partner_pinochle} (partner's pinochle) + {partner_nines} (partner's two 9s) + {all_counters} (all counters) + {last_trick_bonus} (last trick bonus) = {total_score}")
<<<146>>>