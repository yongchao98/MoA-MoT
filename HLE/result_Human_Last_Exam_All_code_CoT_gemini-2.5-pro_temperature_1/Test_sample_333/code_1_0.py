# The horizontal distance from the gun to the highest point of elevation.
I = 500

# Based on the law of conservation of momentum, the second fragment
# will travel a horizontal distance from the point of explosion (the apex)
# that is three times the initial horizontal distance to the apex.
# Let d2 be the distance the second fragment travels from the apex.
d2 = 3 * I

# The maximum safe distance is the total horizontal distance the second
# fragment travels from the gun. This is the sum of the distance to the
# apex and the distance it travels after the explosion.
total_distance = I + d2

print("The problem is solved using the principle of conservation of momentum.")
print(f"The horizontal distance to the apex is I = {I} m.")
print(f"The horizontal distance the second fragment travels after the explosion is 3 * I = {d2} m.")
print("The maximum safe distance is the sum of these two distances.")
print(f"Final calculation: {I} m + {d2} m = {total_distance} m.")