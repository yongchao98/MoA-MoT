import collections

# Step 1: Define the given data and determine the states of the world.
# P(e, y) is the given probability distribution.
prob_dist = {
    (22, 132): 0.4375,
    (22, 44): 0.0625,
    (10, 60): 0.0625,
    (10, 20): 0.4375
}

# The effort levels are the first element of the keys.
efforts = sorted(list(set(e for e, y in prob_dist.keys())), reverse=True)
e_H = efforts[0] # High effort, e_H = 22
e_L = efforts[1] # Low effort, e_L = 10

# Calculate the states 's' from y = s * e.
# We expect two states, s_H and s_L.
states = set()
for (e, y), prob in prob_dist.items():
    s = y / e
    states.add(s)

s_H = max(states) # High state, s_H = 6
s_L = min(states) # Low state, s_L = 2

print(f"From the data, we identify:")
print(f"High effort e_H = {e_H}, Low effort e_L = {e_L}")
print(f"High state s_H = {s_H}, Low state s_L = {s_L}\n")

# Step 2: Map the (e,y) probabilities to (e,s) probabilities
prob_es = {
    (e_H, s_H): prob_dist[(e_H, e_H * s_H)],
    (e_H, s_L): prob_dist[(e_H, e_H * s_L)],
    (e_L, s_H): prob_dist[(e_L, e_L * s_H)],
    (e_L, s_L): prob_dist[(e_L, e_L * s_L)]
}

# Step 3: Calculate conditional expectations E[s|e].
# The effort 'e' serves as a proxy for the unobserved signal 'theta'.

# Marginal probability of high effort P(e_H)
prob_eH = prob_es[(e_H, s_H)] + prob_es[(e_H, s_L)]
# Marginal probability of low effort P(e_L)
prob_eL = prob_es[(e_L, s_H)] + prob_es[(e_L, s_L)]

# Conditional probability P(s|e_H)
prob_sH_given_eH = prob_es[(e_H, s_H)] / prob_eH
prob_sL_given_eH = prob_es[(e_H, s_L)] / prob_eH

# Conditional probability P(s|e_L)
prob_sH_given_eL = prob_es[(e_L, s_H)] / prob_eL
prob_sL_given_eL = prob_es[(e_L, s_L)] / prob_eL

# Conditional expectation E[s|e_H]
E_s_given_eH = s_H * prob_sH_given_eH + s_L * prob_sL_given_eH
# Conditional expectation E[s|e_L]
E_s_given_eL = s_H * prob_sH_given_eL + s_L * prob_sL_given_eL

print("Calculated conditional expectations:")
print(f"E[s | e=e_H] = E[s | e={e_H}] = {E_s_given_eH}")
print(f"E[s | e=e_L] = E[s | e={e_L}] = {E_s_given_eL}\n")

# Step 4: Determine beta from the employee's First Order Condition: e = beta * E[s|e]
# Using the high effort case:
beta_from_eH = e_H / E_s_given_eH

# Using the low effort case:
beta_from_eL = e_L / E_s_given_eL

print("Step 4: Calculate beta (β) from the employee's optimization problem.")
print("The employee's optimal effort rule is: e = β * E[s|θ]")

print("\nFor the high effort case:")
print(f"  {e_H} = β * {E_s_given_eH}")
print(f"  β = {e_H} / {E_s_given_eH} = {beta_from_eH}")

print("\nFor the low effort case:")
print(f"  {e_L} = β * {E_s_given_eL}")
print(f"  β = {e_L} / {E_s_given_eL} = {beta_from_eL}")

# Step 5 & 6: Relate beta to p and find the final answer.
# Since beta values are consistent, we have found the correct beta.
beta = beta_from_eH

print(f"\nThe value of β is consistent across both cases: β = {beta}")
print("\nStep 5 & 6: Find p.")
print("In a principal-agent model with a risk-neutral agent, the firm maximizes profit by setting the piece rate β equal to the price p.")
print("The final equation is: p = β")
print(f"Therefore, the value of p is {beta}.")