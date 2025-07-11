import math

# 1. Define the problem: Find all possible normalized distances r <= 3
#    for hard spheres packed on a 2D plane.

# 2. Model the system using a hexagonal lattice. The squared distance r^2
#    from a central sphere (at the origin) to any other sphere at lattice
#    coordinates (n, m) is given by the formula: r^2 = n^2 + m^2 + n*m.

# 3. We will iterate through integer values for n and m to find all unique
#    distances where r <= 3 (or r^2 <= 9).

# A set is used to automatically store only the unique distances.
distances = set()

# A loop range of -3 to 3 for n and m is sufficient, as any value
# outside this range would result in r^2 > 9.
n_range = range(-3, 4)
m_range = range(-3, 4)

for n in n_range:
    for m in m_range:
        # The distance to the sphere itself is 0, which we ignore.
        if n == 0 and m == 0:
            continue

        # Calculate the squared distance.
        r_squared = n**2 + m**2 + n * m

        # Check if the distance is within the desired range (0 < r <= 3).
        if 0 < r_squared <= 9:
            distance = math.sqrt(r_squared)
            # Add the calculated distance to our set. We round to a high precision
            # to avoid floating-point inaccuracies when checking for uniqueness.
            distances.add(round(distance, 8))

# 4. Sort the unique distances and print them in the required format.
sorted_distances = sorted(list(distances))

print("The set of possible normalized distances r <= 3 is:")
for r in sorted_distances:
    # To show the "equation", we calculate r^2 from r.
    # Rounding ensures we get a clean integer (e.g., 2.999... -> 3).
    r_squared_val = int(round(r**2))
    print("r = sqrt({}) = {:.2f}".format(r_squared_val, r))
