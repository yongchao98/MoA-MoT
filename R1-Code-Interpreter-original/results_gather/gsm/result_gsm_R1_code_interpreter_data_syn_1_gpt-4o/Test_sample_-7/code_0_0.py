# Given dimensions and values
base_length = 6  # inches
base_width = 6   # inches
height = 15      # inches
jelly_bean_volume = 0.15  # cubic inches
packing_efficiency = 0.8
red_jelly_bean_percentage = 0.3
true_red_jelly_beans = 900

# Calculate the volume of the jar
jar_volume = base_length * base_width * height

# Calculate the effective volume for jelly beans
effective_volume = jar_volume * packing_efficiency

# Calculate the total number of jelly beans
total_jelly_beans = effective_volume / jelly_bean_volume

# Estimate the number of red jelly beans
estimated_red_jelly_beans = total_jelly_beans * red_jelly_bean_percentage

# Calculate the difference between the estimated and true number of red jelly beans
difference = estimated_red_jelly_beans - true_red_jelly_beans

# Output the difference
print(difference)