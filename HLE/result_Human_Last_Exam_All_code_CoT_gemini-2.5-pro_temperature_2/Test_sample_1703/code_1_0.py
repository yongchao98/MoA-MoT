import math

def final_answer():
    """
    This function explains and prints the formula for the area of the triangle T(t).
    """

    # Define the constants from the derivation
    # The area is A(t) = C1 + C2 * delta(t)^2
    # C1 = (225 * sqrt(3)) / 4
    # C2 = (3 * sqrt(3)) / 4
    
    C1_num = 225
    C1_den = 4
    C2_num = 3
    C2_den = 4
    sqrt_3_val = math.sqrt(3)

    print("The area of the triangle T(t) as a function of time t is independent of the hexagon's rotation.")
    print("The triangle remains equilateral for all time t.")
    print("Its area, A(t), is given by the following equation:")
    print("")
    
    # Print the equation structure
    print("  A(t) = C1 + C2 * [delta(t)]^2")
    print("")

    # Explain the components
    print("where:")
    print("  - C1 and C2 are constants.")
    print("  - delta(t) is the displacement of each triangle vertex from the midpoint of its respective hexagon side.")
    print("")

    # Output the numbers in the final equation
    print("The final equation with its numerical values is:")
    print(f"  A(t) = ({C1_num} * sqrt(3)) / {C1_den} + ({C2_num} * sqrt(3)) / {C2_den} * [delta(t)]^2")
    print("")
    
    print("In approximate decimal form:")
    print(f"  A(t) approx = {C1_num * sqrt_3_val / C1_den:.4f} + {C2_num * sqrt_3_val / C2_den:.4f} * [delta(t)]^2")
    print("")

    print("The function delta(t) describes the oscillatory motion of the vertices along the hexagon's sides.")
    print("It is a periodic function with a period of 20 seconds. Letting t' = t % 20, it is defined as:")
    print("  delta(t) = t'             , for 0 <= t' <= 5")
    print("  delta(t) = 10 - t'        , for 5 < t' <= 15")
    print("  delta(t) = t' - 20        , for 15 < t' <= 20")

# Execute the function to print the solution
final_answer()