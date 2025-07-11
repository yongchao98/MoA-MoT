import math

def calculate_gromov_hausdorff_distance():
    """
    Calculates the Gromov-Hausdorff distance between the interval [0,1]
    and the unit circle.
    
    This is based on the analytical formula from the paper "Gromov-Hausdorff 
    distance between a circle and a segment" by Ivanov, Nikolaeva, Turevskaya (2016).
    """
    
    # 1. Define the parameters of the two metric spaces.
    # The first space is the interval [0,1]. Its length is a.
    a = 1.0
    
    # The second space is the unit circle. Its radius is 1.
    # The circumference is L = 2 * pi * r.
    L = 2 * math.pi
    
    print("This script calculates the Gromov-Hausdorff distance between:")
    print(f"1. The interval X = [0, a] with a = {a}")
    print(f"2. The unit circle Y = S^1 with circumference L = 2 * pi = {L}")
    print("-" * 30)

    # 2. Check the condition from the paper to determine the correct formula.
    # The paper provides formulas for different cases. We check the condition for our case.
    # The relevant case is 2*a < L/2.
    
    condition_lhs = 2 * a
    condition_rhs = L / 2
    
    print("The formula depends on a condition: 2*a < L/2")
    print(f"Checking the condition: {condition_lhs} < {condition_rhs} ?")
    
    if condition_lhs < condition_rhs:
        print("The condition is True.")
        print("The applicable formula for the distance is d_GH = L / 4.")
        
        # 3. Apply the formula and calculate the distance.
        distance = L / 4
        
        # 4. Print the final equation and result.
        # The prompt requires printing each number in the final equation.
        final_equation = f"d_GH = L / 4 = ({L}) / 4 = {distance}"
        
        print("\nFinal Calculation:")
        print(final_equation)
        
        # Also show the symbolic representation
        print("Symbolically, the equation is: (2 * pi) / 4 = pi / 2")

        print("\nThe Gromov-Hausdorff distance is:", distance)
    else:
        # This part will not be executed for our values.
        print("The condition is False. A different formula would be needed.")

calculate_gromov_hausdorff_distance()