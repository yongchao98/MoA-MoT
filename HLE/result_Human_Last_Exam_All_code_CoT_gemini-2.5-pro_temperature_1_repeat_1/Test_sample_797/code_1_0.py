# Objective: Calculate the total number of turns (right, left, back) for a given FCBD® sequence.

# --- Step 1: Define the turn contributions for each move based on analysis ---

# Swivel Step (no turns)
ss_r, ss_l, ss_b = 0, 0, 0

# Swivel Step Half Turn (16 counts = 2 repetitions of an 8-count move)
# Each repetition shows Right, Left, and Back once.
reps_ssht = 2
ssht_r = 1 * reps_ssht
ssht_l = 1 * reps_ssht
ssht_b = 1 * reps_ssht

# Sunanda (one right turn part, one left turn part; transitional back not counted)
sunanda_r, sunanda_l, sunanda_b = 1, 1, 0

# Balancing Step (no turns)
bs_r, bs_l, bs_b = 0, 0, 0

# Figure 8 (no turns)
f8_r, f8_l, f8_b = 0, 0, 0

# Barrel Turn (full 360-degree spin shows all profiles)
bt_r, bt_l, bt_b = 1, 1, 1

# --- Step 2: Sum the contributions ---

# List of moves that contribute to the turn count
turning_moves_r = [ssht_r, sunanda_r, bt_r]
turning_moves_l = [ssht_l, sunanda_l, bt_l]
turning_moves_b = [ssht_b, sunanda_b, bt_b]

# Calculate totals
total_right = sum(turning_moves_r)
total_left = sum(turning_moves_l)
total_back = sum(turning_moves_b)


# --- Step 3: Print the final equations and result ---

print("The calculation for each turn direction is as follows:")

# Format the equation string by joining the numbers with ' + '
right_equation = ' + '.join(map(str, turning_moves_r))
left_equation = ' + '.join(map(str, turning_moves_l))
back_equation = ' + '.join(map(str, turning_moves_b))

print(f"Right side turns: {right_equation} = {total_right}")
print(f"Left side turns: {left_equation} = {total_left}")
print(f"Back to audience: {back_equation} = {total_back}")
