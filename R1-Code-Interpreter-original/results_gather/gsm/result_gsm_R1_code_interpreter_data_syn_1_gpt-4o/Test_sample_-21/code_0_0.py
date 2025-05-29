# Number of women wearing each type of heel
women_4_inch = 3
women_2_inch = 8449158

# Total height contributed by each group
total_height_4_inch = women_4_inch * 4
total_height_2_inch = women_2_inch * 2

# Total number of women
total_women = women_4_inch + women_2_inch

# Average height of heels
average_height = (total_height_4_inch + total_height_2_inch) / total_women

print(average_height)