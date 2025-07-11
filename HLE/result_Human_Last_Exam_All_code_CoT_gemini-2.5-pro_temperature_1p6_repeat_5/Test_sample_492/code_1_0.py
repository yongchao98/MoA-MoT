import math

def solve_for_c():
    """
    This function derives and calculates the constant c, which represents the time
    of emergence for the giant connected component in the specified random graph model.
    """
    
    # The analysis shows that for large n, the average degree k of the graph at time t
    # can be approximated by the formula k(t) = t^2 / 3.
    
    # The giant component emerges when the average degree k is equal to 1.
    # Therefore, we need to solve the equation c^2 / 3 = 1 for the time c.

    # The equation is: c^2 / 3 = 1
    # We define the numbers in this equation to display them clearly.
    rhs = 1
    denominator = 3
    
    print("The average degree k at time t is given by the formula: k(t) = t^2 / 3")
    print("The giant component emerges when the average degree k becomes 1.")
    print("We solve for the critical time c by setting k(c) = 1.\n")
    
    # Step 1: c^2 / 3 = 1
    print("Step-by-step solution of the final equation:")
    print(f"c^2 / {denominator} = {rhs}")
    
    # Step 2: c^2 = 3
    c_squared = rhs * denominator
    print(f"c^2 = {c_squared}")

    # Step 3: c = sqrt(3)
    c = math.sqrt(c_squared)
    print(f"c = sqrt({c_squared})")
    
    print(f"\nThe exact value of c is sqrt({int(c_squared)}).")
    print(f"The numerical value of c is approximately: {c}")

if __name__ == "__main__":
    solve_for_c()
