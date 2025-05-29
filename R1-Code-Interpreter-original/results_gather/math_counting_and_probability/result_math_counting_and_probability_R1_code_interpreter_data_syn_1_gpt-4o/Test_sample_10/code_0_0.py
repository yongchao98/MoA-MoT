# Number of green and red balls
green_balls = 6
red_balls = 4
total_balls = green_balls + red_balls

# Probability of drawing two green balls
prob_two_green = (green_balls / total_balls) * (green_balls / total_balls)

# Probability of drawing two red balls
prob_two_red = (red_balls / total_balls) * (red_balls / total_balls)

# Total probability of drawing two balls of the same color
prob_same_color = prob_two_green + prob_two_red

print(prob_same_color)