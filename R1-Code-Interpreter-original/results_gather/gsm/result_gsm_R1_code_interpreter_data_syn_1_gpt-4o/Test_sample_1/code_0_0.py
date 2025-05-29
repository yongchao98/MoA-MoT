# Given data
skipping_speed = 3  # in miles per hour
running_speed = skipping_speed * 2  # twice the skipping speed
walking_speed = running_speed / 4  # one-fourth of the running speed

# Time spent on each activity
total_time = 6  # in hours
running_time = total_time / 3  # one-third of the total time
walking_time = total_time * 2 / 3  # two-thirds of the total time

# Calculate distances
distance_running = running_speed * running_time
distance_walking = walking_speed * walking_time

# Total distance
total_distance = distance_running + distance_walking

print(total_distance)