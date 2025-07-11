import sys

def solve_equilibrium():
  """
  Determines if a body is in equilibrium based on force vectors.
  """
  try:
    # Read the number of force vectors
    n = int(sys.stdin.readline())
    
    # Initialize the sum of forces for each coordinate to zero
    sum_x, sum_y, sum_z = 0, 0, 0
    
    # Loop n times to read each vector and add it to the total sum
    for _ in range(n):
      line = sys.stdin.readline()
      # Read the three integer coordinates from the line
      parts = line.split()
      x_i = int(parts[0])
      y_i = int(parts[1])
      z_i = int(parts[2])
      
      # Add the components of the current vector to the sums
      sum_x += x_i
      sum_y += y_i
      sum_z += z_i
      
    # A body is in equilibrium if the sum of all force vectors is the zero vector.
    # The final equation is: sum_x == 0 and sum_y == 0 and sum_z == 0
    # Here are the numbers used in that final equation:
    # print(f"Final vector sum: ({sum_x}, {sum_y}, {sum_z})")

    if sum_x == 0 and sum_y == 0 and sum_z == 0:
      print("YES")
    else:
      print("NO")
      
  except (IOError, ValueError):
    # Assuming the input format is always correct as per the problem statement.
    print("Invalid input.")

solve_equilibrium()