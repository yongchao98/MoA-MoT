# Step 1: Define the points from my meld.
my_meld_run = 15
my_meld_aces = 10

# Step 2: Define the points from my partner's meld.
partner_meld_pinochle = 4
partner_meld_dix = 2

# Step 3: Define the points from tricks.
# In 48-card pinochle, there are 24 counter cards (A, 10, K), each worth 1 point.
# Having all 8 aces and playing perfectly guarantees winning all tricks.
trick_points = 24

# Step 4: Calculate the total points.
total_points = my_meld_run + my_meld_aces + partner_meld_pinochle + partner_meld_dix + trick_points

# Step 5: Print the breakdown and the final answer.
print("Calculating the total score for the hand:")
print("My meld points:")
print(f"  - Run: {my_meld_run} points")
print(f"  - 8 Aces: {my_meld_aces} points")
print("Partner's meld points:")
print(f"  - Pinochle: {partner_meld_pinochle} points")
print(f"  - Two 9s of trump: {partner_meld_dix} points")
print("Trick points:")
print(f"  - Points from all 12 tricks: {trick_points} points")
print("\nThe final score is the sum of all these points.")
print(f"Final Equation: {my_meld_run} + {my_meld_aces} + {partner_meld_pinochle} + {partner_meld_dix} + {trick_points} = {total_points}")
print(f"Total points earned: {total_points}")
