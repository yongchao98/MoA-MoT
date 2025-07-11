import sympy

# Let's define the probabilities for the edges.
# The probability that a vertical edge is deleted is 1/2.
p_del_vert = sympy.Rational(1, 2)
p_exist_vert = 1 - p_del_vert

# The probability that an upper horizontal edge is deleted is 1/3.
p_del_horz_up = sympy.Rational(1, 3)
p_exist_horz_up = 1 - p_del_horz_up

print("Step 1: Define transition rates between levels as c -> oo")
# R_01: Rate from level 0 to 1.
# A walker at (n,0) moves to (n,1) with probability proportional to e^0 = 1,
# while moving right to (n+1,0) is proportional to e^c.
# The transition probability P((n,0)->(n,1)) is roughly e^(-c).
# R_01 is the average of this probability. The leading term for large c is:
# R_01 ~ (1/2) * e^(-c). So, lim_{c->oo} R_01 = 0.

# R_10: Rate from level 1 to 0.
# A walker at (n,1) will be forced to consider a vertical move if the
# upper horizontal edge to (n+1,1) is missing (prob 1/3).
# If the vertical edge to (n,0) exists (prob 1/2), this move has weight 1,
# which is infinitely larger than moving left (weight e^(-c)).
# We need to average over local edge configurations for a node (n,1).
# A walker can only exist on the infinite component. A node (n,1) is isolated
# if its edges to (n+1,1), (n-1,1), and (n,0) are all missing.
prob_isolated_node_at_1 = p_del_horz_up * p_del_horz_up * p_del_vert
prob_on_infinite_component = 1 - prob_isolated_node_at_1

# The transition 1->0 is taken with probability near 1 if the edge (n,1)->(n+1,1)
# is missing AND the edge (n,1)->(n,0) exists.
prob_force_down = p_del_horz_up * p_exist_vert
# We need to normalize by the probability of not being on an isolated node.
R_10_limit = prob_force_down / prob_on_infinite_component

print(f"The limit of the transition rate R_10 as c -> infinity is: {prob_force_down} / {prob_on_infinite_component} = {R_10_limit}")

# As R_01 -> 0 and R_10 -> a positive constant, the ratio pi_1/pi_0 = R_01/R_10 -> 0.
# This means pi_0 -> 1 and pi_1 -> 0.
print("\nStep 2: Determine the stationary distribution as c -> oo")
print("The walker spends almost all its time on the lower level (level 0).")
pi_0_limit = 1
pi_1_limit = 0
print(f"lim pi_0 = {pi_0_limit}")
print(f"lim pi_1 = {pi_1_limit}")


print("\nStep 3: Calculate the speed on each level as c -> oo")
# v_0: Speed on level 0.
# The lower horizontal edge ((n,0),(n+1,0)) always exists. The walker at (n,0)
# sees a rightward move with weight e^c, which is infinitely preferred.
# So the walker always moves right. One step horizontally in one unit of time.
v_0_limit = 1
print(f"The limiting speed on level 0, v_0, is: {v_0_limit}")

# v_1: Speed on level 1.
# This is the expected horizontal displacement per step from a node on level 1.
# A step is rightward (+1) if the edge ((n,1),(n+1,1)) exists (prob 2/3).
# A step is leftward (-1) if the rightward edge is missing (prob 1/3), the downward
# edge is missing (prob 1/2), and the leftward edge exists (prob 2/3).
prob_move_right = p_exist_horz_up
prob_move_left = p_del_horz_up * p_del_vert * p_exist_horz_up

# Average displacement, normalized by being on the infinite component.
v_1_limit = (prob_move_right * 1 + prob_move_left * (-1)) / prob_on_infinite_component
print(f"The limiting speed on level 1, v_1, is: ({prob_move_right}*1 + {prob_move_left}*(-1)) / {prob_on_infinite_component} = {v_1_limit}")

print("\nStep 4: Combine results to find the final asymptotic speed")
# The final speed is v = pi_0*v_0 + pi_1*v_1
# We take the limit as c -> oo
final_speed = pi_0_limit * v_0_limit + pi_1_limit * v_1_limit
print(f"The limit of the asymptotic speed v(c) is lim(pi_0*v_0 + pi_1*v_1)")
print(f"= (lim pi_0)(lim v_0) + (lim pi_1)(lim v_1)")
print(f"= {pi_0_limit} * {v_0_limit} + {pi_1_limit} * {v_1_limit}")
print(f"= {final_speed}")