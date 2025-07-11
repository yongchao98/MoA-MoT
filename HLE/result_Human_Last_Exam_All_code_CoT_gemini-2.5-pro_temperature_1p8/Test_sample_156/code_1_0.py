import math

def solve_alien_colony():
    """
    This function determines the maximal size of the alien colony.
    It follows the logical steps outlined above.
    """
    
    # Step 1 & 2: Determine the minimum initial perimeter for an 8-square colony.
    # The most compact shapes (e.g., 4x2 rectangle) have a perimeter of 12.
    initial_area = 8
    min_initial_perimeter = 12

    print(f"The colony starts with {initial_area} squares.")
    print("The expansion rule means the perimeter of the colony can never increase.")
    print(f"The minimum possible perimeter for the initial {initial_area} squares is {min_initial_perimeter}.")
    print("The final stable shape must be a rectangle.")
    print(f"Therefore, we seek the largest rectangle with a perimeter of {min_initial_perimeter} or less.\n")
    
    # Step 4: Find the rectangle with max area given the perimeter constraint.
    # A rectangle's perimeter is 2*(width + height).
    # 2*(w+h) <= 12  =>  w+h <= 6
    max_area = 0
    best_dims = (0, 0)

    # We iterate through possible rectangle dimensions (w,h) on an 8x8 board.
    print("Checking possible stable rectangle dimensions:")
    for w in range(1, 9):
      for h in range(w, 9): # Iterate h from w to avoid duplicate (w,h) and (h,w) pairs
          perimeter = 2 * (w + h)
          if perimeter <= min_initial_perimeter:
              area = w * h
              print(f" - {w}x{h} rectangle: Area = {area}, Perimeter = {perimeter}")
              if area > max_area:
                  max_area = area
                  best_dims = (w, h)
                  
    print(f"\nThe largest possible stable area is {max_area}, from a {best_dims[0]}x{best_dims[1]} rectangle.\n")
    
    # Step 5: Confirming that K is achievable.
    print("This maximum size is achievable.")
    print("Start with an 8-square colony shaped like a 3x3 square with one corner missing.")
    print("This arrangement has a perimeter of 12 and can contain the required d5 and e5 squares.")
    print("The vacant corner square has 2 captured neighbors, so it gets captured.")
    print("The colony grows into a stable 3x3 square, reaching the maximal size.")
    
    # Final Result
    final_w, final_h = best_dims
    print("\nThe equation for the final maximal size K is:")
    print(f"{final_w} * {final_h} = {max_area}")

solve_alien_colony()
<<<9>>>