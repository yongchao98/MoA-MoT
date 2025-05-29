from math import comb

# Total number of men and women
men = 3
women = 4

# Scenario 1: Group of 3 has 1 man and 2 women
ways_1_man_2_women = comb(men, 1) * comb(women, 2)

# Remaining: 2 men and 2 women
# Form two groups of 2, each with 1 man and 1 woman
ways_2_groups_1_man_1_woman = comb(2, 1) * comb(2, 1)

# Total ways for scenario 1
total_ways_scenario_1 = ways_1_man_2_women * ways_2_groups_1_man_1_woman

# Scenario 2: Group of 3 has 2 men and 1 woman
ways_2_men_1_woman = comb(men, 2) * comb(women, 1)

# Remaining: 1 man and 3 women
# Form two groups of 2, each with 1 man and 1 woman
ways_2_groups_1_man_1_woman_scenario_2 = comb(1, 1) * comb(3, 1)

# Total ways for scenario 2
total_ways_scenario_2 = ways_2_men_1_woman * ways_2_groups_1_man_1_woman_scenario_2

# Total ways
total_ways = total_ways_scenario_1 + total_ways_scenario_2

print(total_ways)