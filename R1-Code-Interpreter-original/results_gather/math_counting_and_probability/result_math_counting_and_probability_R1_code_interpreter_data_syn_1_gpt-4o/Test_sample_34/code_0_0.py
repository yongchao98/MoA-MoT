# Probabilities of selecting each container
p_container = 1/3

# Probabilities of selecting a green ball from each container
p_green_I = 4/12
p_green_II = 4/6
p_green_III = 4/6

# Total probability of selecting a green ball
p_green = (p_green_I * p_container) + (p_green_II * p_container) + (p_green_III * p_container)

print(p_green)