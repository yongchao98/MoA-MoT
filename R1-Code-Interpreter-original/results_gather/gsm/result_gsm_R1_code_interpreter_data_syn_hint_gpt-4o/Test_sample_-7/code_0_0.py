# Dimensions of the jar
base_length = 6  # inches
base_width = 6   # inches
height = 15      # inches

# Volume of the jar
volume_jar = base_length * base_width * height

# Packing efficiency
packing_efficiency = 0.8

# Effective volume for jelly beans
effective_volume = volume_jar * packing_efficiency

# Volume of a single jelly bean
volume_jelly_bean = 0.15  # cubic inches

# Total number of jelly beans
total_jelly_beans = effective_volume / volume_jelly_bean

# Percentage of red jelly beans
percentage_red = 0.3

# Estimated number of red jelly beans
estimated_red_jelly_beans = total_jelly_beans * percentage_red

# True number of red jelly beans
true_red_jelly_beans = 900

# Difference between estimate and true number
difference = estimated_red_jelly_beans - true_red_jelly_beans

print(difference)