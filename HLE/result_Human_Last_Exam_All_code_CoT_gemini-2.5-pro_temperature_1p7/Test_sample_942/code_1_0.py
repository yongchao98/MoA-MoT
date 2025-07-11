import numpy as np

def l_infinity_norm(v):
  """Calculates the L-infinity norm of a vector."""
  return np.max(np.abs(v))

# Define the endpoints of the line segment
p = np.array([1.0, 1.0])
q = np.array([-1.0, 1.0])

# Verify that p and q are on the unit sphere
norm_p = l_infinity_norm(p)
norm_q = l_infinity_norm(q)

print(f"Point p = {p}")
print(f"Point q = {q}")
print("-" * 30)
print("Verifying points are on the L-infinity unit sphere:")
print(f"||p|| = max(|{p[0]}|, |{p[1]}|) = {norm_p}")
print(f"||q|| = max(|{q[0]}|, |{q[1]}|) = {norm_q}")
print("-" * 30)

# Verify that the line segment lies on the unit sphere by checking a few points
print("Verifying segment is on the unit sphere (checking t in [0,1]):")
for t in np.linspace(0, 1, 5):
    point_on_segment = (1 - t) * p + t * q
    norm_of_point = l_infinity_norm(point_on_segment)
    print(f"t={t:.2f}: point={point_on_segment}, norm = {norm_of_point:.1f}")
print("-" * 30)

# Calculate the length of the segment
p_minus_q = p - q
length = l_infinity_norm(p_minus_q)

print("Calculating the length of the segment ||p - q||:")
print(f"p - q = {p} - {q} = {p_minus_q}")
print(f"Length = ||p - q|| = max(|{p_minus_q[0]}|, |{p_minus_q[1]}|) = {length}")
