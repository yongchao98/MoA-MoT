# Total number of chips
total_chips = 12

# Probability of drawing two blue chips
prob_blue_blue = (5/total_chips) * (5/total_chips)

# Probability of drawing two red chips
prob_red_red = (4/total_chips) * (4/total_chips)

# Probability of drawing two yellow chips
prob_yellow_yellow = (3/total_chips) * (3/total_chips)

# Probability of drawing two chips of the same color
prob_same_color = prob_blue_blue + prob_red_red + prob_yellow_yellow

# Probability of drawing two chips of different colors
prob_different_color = 1 - prob_same_color

# Expressing the probability as a common fraction
from fractions import Fraction
prob_different_color_fraction = Fraction(prob_different_color).limit_denominator()

print(prob_different_color_fraction)