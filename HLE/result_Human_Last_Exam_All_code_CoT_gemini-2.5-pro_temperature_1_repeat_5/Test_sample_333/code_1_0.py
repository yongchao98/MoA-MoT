# Input variable:
# I is the horizontal distance from the gun to the highest point of elevation.
I = 500

# The projectile explodes at its highest point, I = 500 m from the gun.
# The explosion splits it into two equal parts.
# To maximize the range of the second fragment, the first fragment must be
# propelled backward to land at the starting point (the gun).
# This means its horizontal velocity is reversed.

# By conservation of momentum, if the first fragment's velocity becomes -vx,
# the second fragment's velocity must become 3*vx to conserve the total momentum.
# (M*vx = (M/2)*(-vx) + (M/2)*(3*vx))

# The time to fall from the apex is the same as the time to rise.
# The original distance to the apex is I = vx * t_rise.
# The additional distance the second fragment travels after the explosion is d2 = (3*vx) * t_fall.
# Since t_rise = t_fall, d2 = 3 * I.
distance_traveled_by_second_fragment = 3 * I

# The total distance from the gun is the distance to the apex plus the additional distance
# traveled by the second fragment.
total_distance = I + distance_traveled_by_second_fragment

# Print the final equation with the numbers.
print(f"The maximum distance from the gun is the sum of the distance to the apex ({I} m) and the additional distance the second fragment travels ({int(distance_traveled_by_second_fragment)} m).")
print("Final Equation:")
print(f"{I} + 3 * {I} = {int(total_distance)}")
