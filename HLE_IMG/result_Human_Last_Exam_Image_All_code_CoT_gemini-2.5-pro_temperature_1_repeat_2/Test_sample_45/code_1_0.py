def find_rotational_symmetry():
  """
  This function explains the process of finding the rotational symmetry of the provided tiling.
  """
  print("1. Rotational symmetry of order 'n' means the pattern looks the same after being rotated by 360/n degrees.")
  print("2. We look for points in the tiling that can act as centers of rotation.")
  print("3. By observing the pattern, we can identify centers of 6-fold symmetry at the center of each blue hexagon.")
  print("4. We can also identify centers of 4-fold symmetry at the center of the yellow squares in the star-like motifs.")
  print("5. The rotational symmetry of the entire tiling is defined by the highest order of symmetry found.")
  
  highest_order_symmetry = 6
  full_circle = 360
  rotation_angle = full_circle / highest_order_symmetry
  
  print(f"\nThe highest order of rotational symmetry found is {highest_order_symmetry}.")
  print("This corresponds to a rotation angle that leaves the pattern unchanged.")
  print(f"The calculation for this angle is: {full_circle} / {highest_order_symmetry} = {int(rotation_angle)} degrees.")
  print(f"Therefore, the rotational symmetry of the tiling is {highest_order_symmetry}.")

find_rotational_symmetry()