import math

# Given values
c_hb = 0.9375  # cos(<h_p, b_p>)
c_hz = 0.9    # cos(<h_p, z_p>)

# The relationship between the cosines of the angles for coplanar vectors
# can be expressed in two ways, depending on their geometric arrangement.
# A common physical assumption in such models is that the resulting eigenvector h_p
# lies on the geodesic connecting the underlying signal directions b_p and z_p.
# This corresponds to the angular relation: angle(b,z) = angle(b,h) + angle(h,z).
# The cosine of this sum of angles is given by the formula:
# cos(A + B) = cos(A)cos(B) - sin(A)sin(B)

# Calculate the sines
# s_hb = sin(<h_p, b_p>)
s_hb = math.sqrt(1 - c_hb**2)

# s_hz = sin(<h_p, z_p>)
s_hz = math.sqrt(1 - c_hz**2)

# Calculate the final value c_bz = lim <b_p, z_p>
c_bz = c_hb * c_hz - s_hb * s_hz

# Print the components of the calculation
print(f"Let c_hb = lim <h_p, b_p> and c_hz = lim <h_p, z_p>.")
print(f"Let c_bz = lim <b_p, z_p>, the value we want to find.")
print(f"Given c_hb = {c_hb}")
print(f"Given c_hz = {c_hz}")
print(f"We calculate s_hb = sqrt(1 - c_hb^2) = {s_hb}")
print(f"We calculate s_hz = sqrt(1 - c_hz^2) = {s_hz}")
print("\nAssuming the vectors are coplanar and h_p lies on the geodesic between b_p and z_p, we use the identity:")
print("c_bz = c_hb * c_hz - s_hb * s_hz")
print("\nSubstituting the values into the equation:")
# The final expression includes all the numbers used
print(f"c_bz = {c_hb} * {c_hz} - {s_hb} * {s_hz}")

# Print the final result
print(f"\nThe value of lim <b_p, z_p> is: {c_bz}")

print("\nFinal Answer:")
print(f"<<<{c_bz}>>>")