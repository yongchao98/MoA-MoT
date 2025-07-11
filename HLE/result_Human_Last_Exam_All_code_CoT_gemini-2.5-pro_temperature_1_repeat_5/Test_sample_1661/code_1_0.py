def check_equilibrium():
  """
  Determines if a body is in equilibrium based on force vectors.
  """
  try:
    # 1. Read the number of vectors.
    n = int(input())
    
    # 2. Initialize sums for each coordinate to zero.
    sum_x = 0
    sum_y = 0
    sum_z = 0
    
    # 3. Loop n times to read each vector and add its components to the sums.
    for _ in range(n):
      line = input().split()
      x = int(line[0])
      y = int(line[1])
      z = int(line[2])
      
      sum_x += x
      sum_y += y
      sum_z += z
      
    # 4. Check if the total force vector is the zero vector.
    if sum_x == 0 and sum_y == 0 and sum_z == 0:
      print("YES")
    else:
      print("NO")

  except (ValueError, IndexError):
    print("Error: Invalid input format. Please provide integers as specified.")

check_equilibrium()