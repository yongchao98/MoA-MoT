# Number of plants
total_plants = 20
plants_half_cup = 4
plants_one_cup = 6
plants_quarter_cup = total_plants - plants_half_cup - plants_one_cup

# Water needed
water_half_cup = plants_half_cup * 0.5
water_one_cup = plants_one_cup * 1
water_quarter_cup = plants_quarter_cup * 0.25

# Total water needed
total_water = water_half_cup + water_one_cup + water_quarter_cup
print(total_water)