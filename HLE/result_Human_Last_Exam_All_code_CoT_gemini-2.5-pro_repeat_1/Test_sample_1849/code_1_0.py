import numpy as np

# This script demonstrates one step of the Value Iteration algorithm.
# It shows how the new value function V_{k+1} is calculated from V_k.

# MDP Parameters
# Two states: S0, S1
# V_k: The value function at iteration k
V_k = {'S0': 10.0, 'S1': -1.0}

# Rewards R(s, a)
# For S0: R(S0, 'a') = 5, R(S0, 'b') = 10
# For S1: R(S1, 'c') = -1
rewards = {
    'S0': {'a': 5, 'b': 10},
    'S1': {'c': -1}
}

# Transition Probabilities P(s' | s, a)
# From S0, action 'a': 50% to S0, 50% to S1
# From S0, action 'b': 100% to S1
# From S1, action 'c': 100% to S1
transitions = {
    'S0': {
        'a': {'S0': 0.5, 'S1': 0.5},
        'b': {'S0': 0.0, 'S1': 1.0}
    },
    'S1': {
        'c': {'S0': 0.0, 'S1': 1.0}
    }
}

# Discount Factor
gamma = 0.9

# --- Calculate V_{k+1}(S1) ---
# S1 has only one action 'c'
r_s1_c = rewards['S1']['c']
p_s1_s1_c = transitions['S1']['c']['S1']
v_k_s1 = V_k['S1']
v_kplus1_s1 = r_s1_c + gamma * (p_s1_s1_c * v_k_s1)

print("Calculating V_{k+1} for State S1:")
print(f"V(S1) = R(S1,'c') + γ * P(S1|S1,'c') * V_k(S1)")
print(f"V(S1) = {r_s1_c} + {gamma} * {p_s1_s1_c} * {v_k_s1}")
print(f"V(S1) = {v_kplus1_s1:.2f}\n")


# --- Calculate V_{k+1}(S0) ---
# S0 has two actions 'a' and 'b'. We must calculate the value for each and take the max.

# Action 'a'
r_s0_a = rewards['S0']['a']
p_s0_s0_a = transitions['S0']['a']['S0']
p_s0_s1_a = transitions['S0']['a']['S1']
v_k_s0 = V_k['S0']
q_s0_a = r_s0_a + gamma * (p_s0_s0_a * v_k_s0 + p_s0_s1_a * v_k_s1)

print("Calculating V_{k+1} for State S0:")
print("Value of action 'a':")
print(f"Q(S0,'a') = R(S0,'a') + γ * (P(S0|S0,'a')*V_k(S0) + P(S1|S0,'a')*V_k(S1))")
print(f"Q(S0,'a') = {r_s0_a} + {gamma} * ({p_s0_s0_a} * {v_k_s0} + {p_s0_s1_a} * {v_k_s1})")
print(f"Q(S0,'a') = {q_s0_a:.2f}\n")

# Action 'b'
r_s0_b = rewards['S0']['b']
p_s0_s1_b = transitions['S0']['b']['S1']
q_s0_b = r_s0_b + gamma * (p_s0_s1_b * v_k_s1)

print("Value of action 'b':")
print(f"Q(S0,'b') = R(S0,'b') + γ * P(S1|S0,'b') * V_k(S1)")
print(f"Q(S0,'b') = {r_s0_b} + {gamma} * {p_s0_s1_b} * {v_k_s1}")
print(f"Q(S0,'b') = {q_s0_b:.2f}\n")

# V_{k+1}(S0) is the max over the actions
v_kplus1_s0 = max(q_s0_a, q_s0_b)

print("V_{k+1}(S0) is the maximum of the action values:")
print(f"V(S0) = max({q_s0_a:.2f}, {q_s0_b:.2f})")
print(f"V(S0) = {v_kplus1_s0:.2f}")