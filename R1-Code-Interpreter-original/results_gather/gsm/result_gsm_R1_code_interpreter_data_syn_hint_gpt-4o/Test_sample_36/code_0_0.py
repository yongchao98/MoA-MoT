# Number of plants
total_plants = 20
half_cup_plants = 4
one_cup_plants = 6  # Assumed number due to error in the problem statement
quarter_cup_plants = total_plants - half_cup_plants - one_cup_plants

# Water needed
water_for_half_cup_plants = half_cup_plants * 0.5
water_for_one_cup_plants = one_cup_plants * 1
water_for_quarter_cup_plants = quarter_cup_plants * 0.25

# Total water needed
total_water = water_for_half_cup_plants + water_for_one_cup_plants + water_for_quarter_cup_plants

print(total_water)