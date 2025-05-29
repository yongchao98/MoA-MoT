# Probabilities of selecting a green ball from each container
p_green_A = 3/5
p_green_B = 2/5
p_green_C = 2/5

# Probability of selecting each container
p_container = 1/3

# Total probability of selecting a green ball
p_green = p_container * p_green_A + p_container * p_green_B + p_container * p_green_C

# Output the result
print(p_green)