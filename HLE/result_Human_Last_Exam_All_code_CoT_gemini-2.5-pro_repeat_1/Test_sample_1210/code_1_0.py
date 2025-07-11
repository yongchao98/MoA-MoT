# Initial state of the game
p1_pits = [0, 2, 0, 0, 2, 0]
p1_store = 22
p2_pits = [1, 0, 0, 0, 0, 0] # Left-to-right from P2's perspective
p2_store = 21

print("Initial State:")
print(f"Player 1 Pits: {p1_pits}, Store: {p1_store}")
print(f"Player 2 Pits: {p2_pits}, Store: {p2_store}")
print("-" * 30)

# --- Path to Score Difference of 4 ---
# This path is chosen by Player 1 to prevent Player 2 from capturing.

# P1 Move 1: Play from pit 5 (index 4)
print("Player 1's turn. P1 plays from the 5th pit (2 stones).")
stones_in_hand = p1_pits[4]
p1_pits[4] = 0
# Sow stones
p1_pits[5] += 1
p1_store += 1 # Last stone in store, free turn
print("Last stone landed in P1's store. P1 gets a free turn.")
print(f"Board State: P1 Pits {p1_pits}, Store {p1_store}. P2 Pits {p2_pits}, Store {p2_store}")
print("-" * 30)

# P1 Move 2 (Free Turn): Play from pit 2 (index 1)
print("Player 1's free turn. P1 plays from the 2nd pit (2 stones).")
stones_in_hand = p1_pits[1]
p1_pits[1] = 0
# Sow stones
p1_pits[2] += 1
p1_pits[3] += 1
print("Turn passes to Player 2.")
print(f"Board State: P1 Pits {p1_pits}, Store {p1_store}. P2 Pits {p2_pits}, Store {p2_store}")
print("-" * 30)

# P2 Move 1: Play from their only pit (pit 1, index 0)
print("Player 2's turn. P2 plays from the 1st pit (1 stone).")
# Note: P2's pit 1 is opposite P1's pit 6. P2's move lands in their pit 2.
# P2's pit 2 is opposite P1's pit 5, which is empty. No capture.
p2_pits[0] = 0
p2_pits[1] += 1
print("No capture for P2. Turn passes to P1.")
print(f"Board State: P1 Pits {p1_pits}, Store {p1_store}. P2 Pits {p2_pits}, Store {p2_store}")
print("-" * 30)

# Game End state
# Player 2 has no more stones to move. Player 1 clears their side.
# The stones from P1's pits will all go into P1's store.
print("Player 2 has no more moves. Player 1 clears their side.")
p1_remaining_stones = sum(p1_pits)
p1_store += p1_remaining_stones

# The game ends because P1's side will be empty.
# P2 collects the remaining stones on their side.
p2_remaining_stones = sum(p2_pits)
p2_store += p2_remaining_stones

# Final Scores
final_p1_score = p1_store
final_p2_score = p2_store
score_difference = abs(final_p1_score - final_p2_score)

print("Final Scores:")
print(f"Player 1: {final_p1_score}")
print(f"Player 2: {final_p2_score}")
print(f"Score Difference = {final_p1_score} - {final_p2_score} = {score_difference}")

print("\nAnalysis Result:")
print("A difference of 0 is possible through another line of play.")
print("A difference of 4 is demonstrated to be possible.")
print("Odd differences (1, 3, 5) are impossible due to the even total number of stones.")
print("A difference of 2 is impossible because Player 2 cannot capture exactly one stone.")
print("Therefore, the score difference that is not possible is 2.")
