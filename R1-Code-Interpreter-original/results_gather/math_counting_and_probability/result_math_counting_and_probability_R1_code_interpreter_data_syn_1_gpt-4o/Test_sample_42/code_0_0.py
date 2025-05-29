# Probabilities
p_Alice_keeps = 1/2
p_Alice_tosses_to_Bob = 1/2
p_Bob_tosses_to_Alice = 2/5

# Scenario 1: Alice keeps the ball both turns
prob_scenario_1 = p_Alice_keeps * p_Alice_keeps

# Scenario 2: Alice tosses to Bob, Bob tosses back to Alice
prob_scenario_2 = p_Alice_tosses_to_Bob * p_Bob_tosses_to_Alice

# Total probability that Alice has the ball after two turns
total_probability = prob_scenario_1 + prob_scenario_2

print(total_probability)