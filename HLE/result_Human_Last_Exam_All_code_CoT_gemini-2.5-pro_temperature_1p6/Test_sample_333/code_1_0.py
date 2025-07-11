import sys

def solve_projectile_explosion():
    """
    Calculates the landing position of the second fragment of an exploded projectile.
    """
    # Given values from the problem description
    I = 500  # m, horizontal distance to the highest point
    
    # The first fragment fell near the gun, so its final landing position is 0.
    x1_final = 0  # m
    
    # --- Step 1: Determine the landing position of the Center of Mass (CM) ---
    # The explosion is an internal force and does not affect the trajectory of the CM.
    # The CM will follow the original path of the projectile.
    # For a symmetrical trajectory, the total range is twice the distance to the highest point.
    x_cm_final = 2 * I
    
    # --- Step 2: Use the Center of Mass formula to find the second fragment's position ---
    # The formula for the CM of two equal masses is: x_cm = (x1 + x2) / 2
    # We can rearrange this to solve for the landing position of the second fragment (x2):
    # x2 = 2 * x_cm - x1
    x2_final = 2 * x_cm_final - x1_final
    
    # --- Step 3: Print the explanation and results ---
    print("This problem can be solved by considering the motion of the center of mass (CM) of the two fragments.")
    print("The explosion, being an internal force, does not alter the path of the CM.")
    print("\n1. The projectile reaches its highest point at a horizontal distance I = {} m.".format(I))
    print("2. The CM of the fragments will continue on this path, landing at a total range of 2 * I = {} m.".format(x_cm_final))
    print("3. One fragment lands near the gun, so its final position is x1 = {} m.".format(x1_final))
    print("\nFor two fragments of equal mass, the landing position of the CM (x_cm) is the average of their individual landing positions (x1, x2):")
    print("Equation: x_cm = (x1 + x2) / 2")
    
    print("\nTo find the landing position of the second fragment (x2), we solve for it:")
    print("x2 = 2 * x_cm - x1")
    print("Plugging in the values, we get:")
    print("x2 = 2 * {} - {}".format(x_cm_final, x1_final))
    
    print("\nThe final equation with all numbers substituted is:")
    print("{:.0f} = ({:.0f} + {:.0f}) / 2".format(x_cm_final, x1_final, x2_final))
    
    print("\nTherefore, the maximum distance from the gun you can safely be is {} meters.".format(int(x2_final)))
    
    # Returning the final answer for the wrapper
    return int(x2_final)

# Execute the function and capture the answer
final_answer = solve_projectile_explosion()

# The final line of the output must be the answer in the specified format
sys.stdout.write(f"\n<<<{final_answer}>>>\n")
