# Given data
base_length = 6  # inches
base_width = 6   # inches
height = 15      # inches
jelly_bean_volume = 0.15  # cubic inches
packing_efficiency = 0.80
red_percentage = 0.30
actual_red_jelly_beans = 900

# Calculate the volume of the jar
jar_volume = base_length * base_width * height

# Calculate the effective volume for jelly beans
effective_volume = jar_volume * packing_efficiency

# Calculate the total number of jelly beans
total_jelly_beans = effective_volume / jelly_bean_volume

# Calculate the estimated number of red jelly beans
estimated_red_jelly_beans = total_jelly_beans * red_percentage

# Calculate how far off Tim's guess was
difference = estimated_red_jelly_beans - actual_red_jelly_beans

print(difference)