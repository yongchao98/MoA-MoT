import math

def calculate_giant_component_time():
    """
    This function explains and calculates the constant 'c' for the emergence
    of the giant component in the described random graph model.
    """
    print("This program calculates the constant 'c', the time of emergence for the giant component.")
    
    print("\nStep 1: State the condition for the emergence of the giant component.")
    print("The giant component emerges when the average vertex degree is 1.")
    print("The condition is: N * p = 1")
    print("  - N is the number of vertices.")
    print("  - p is the probability of an edge between any two vertices.\n")
    
    print("Step 2: Define N and p at the critical time 'c'.")
    print("At time c, the number of vertices is N ≈ n*c.")
    print("The edge probability, for large n, is p ≈ c / (3*n).\n")

    print("Step 3: Substitute N and p into the condition and solve for c.")
    print("The equation becomes: (n * c) * (c / (3 * n)) = 1")
    print("The 'n' terms cancel out, simplifying the equation:\n")
    
    # Define the numbers in the final equation
    lhs_numerator = 1
    lhs_denominator = 3
    rhs = 1
    
    print(f"   c^2")
    print(f"-------- = {rhs}")
    print(f"    {lhs_denominator}\n")
    
    print("Now, we solve for c^2:")
    # Calculate c_squared
    c_squared = rhs * lhs_denominator
    print(f"c^2 = {rhs} * {lhs_denominator}")
    print(f"c^2 = {c_squared}\n")

    # Calculate c
    c = math.sqrt(c_squared)
    print("Finally, we take the square root to find c:")
    print(f"c = sqrt({c_squared})")
    
    print("\n----------------------------------------------------")
    print(f"The exact value of c is sqrt(3).")
    print(f"The numerical value of c is: {c}")
    print("----------------------------------------------------")

# Execute the function
calculate_giant_component_time()