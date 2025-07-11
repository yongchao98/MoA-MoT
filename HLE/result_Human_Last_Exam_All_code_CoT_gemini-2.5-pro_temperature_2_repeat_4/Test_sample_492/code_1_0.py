import math

def find_critical_time_c():
    """
    Calculates and explains the derivation of 'c', the time of emergence 
    of the giant connected component in the specified random graph model.
    """
    print("Here is the step-by-step derivation for the value of 'c':")
    print("-" * 50)
    print("1. Expected number of vertices at time t, E[V(t)]:")
    print("   Vertices arrive at rate n, so E[V(t)] = n * t.")
    print("\n2. Expected number of edges at time t, E[E(t)]:")
    print("   The rate of edge addition at time s is (1/n) * E[V(s) choose 2] = (1/n) * (n*s)^2 / 2 = n*s^2 / 2.")
    print("   Integrating this rate from 0 to t gives E[E(t)] = integral(n*s^2 / 2) ds = n * t^3 / 6.")
    print("\n3. Average degree d(t) at time t:")
    print("   d(t) = 2 * E[E(t)] / E[V(t)]")
    print("   d(t) = 2 * (n * t^3 / 6) / (n * t) = t^2 / 3.")
    print("\n4. Critical condition for giant component emergence:")
    print("   The giant component emerges when the average degree d(c) = 1.")

    # The final equation is c^2 / 3 = 1
    # We define the numbers in this equation to print them explicitly.
    c_power = 2
    denominator = 3
    rhs = 1  # Right Hand Side

    print("\nThis gives the final equation to solve for the critical time 'c':")
    print(f"c^{c_power} / {denominator} = {rhs}")

    # Solve the equation: c^2 = 3 * 1
    solution_squared = denominator * rhs

    print("\n5. Solving for 'c':")
    print(f"   c^2 = {solution_squared}")
    
    # c must be positive since it represents time.
    final_c_value = math.sqrt(solution_squared)
    final_c_exact = f"sqrt({solution_squared})"
    
    print(f"   c = {final_c_exact}")
    print("\n----------------------------------------------------")
    print(f"The exact value of c is {final_c_exact}.")
    print(f"The approximate numerical value of c is {final_c_value}")
    print("----------------------------------------------------")

# Execute the function to perform the calculation and print the results.
find_critical_time_c()