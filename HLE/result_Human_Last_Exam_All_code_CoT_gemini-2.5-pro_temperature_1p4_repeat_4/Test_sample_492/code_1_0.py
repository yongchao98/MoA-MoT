import math

def solve_for_c():
    """
    This function follows a mean-field approach to find the time 'c'
    of the emergence of the giant component in the described graph model.
    It prints the step-by-step derivation.
    """

    print("Derivation of the time 'c' for the giant component emergence:\n")

    # Step 1: Model the number of vertices V at time c.
    # Vertices arrive at a stochastic rate of n. In the large n limit, we can model this
    # deterministically. The rate of change of vertices is dv/dt = n.
    # Integrating from 0 to c gives the number of vertices V.
    print("Step 1: Determine the number of vertices (V) at the critical time 'c'.")
    print("Vertices arrive at a rate 'n'. So, after time 'c', the number of vertices is:")
    print("V = n * c")
    print("-" * 50)

    # Step 2: Model the number of edges E at time c.
    # Edges form between any pair of existing vertices at a rate of 1/n.
    # At time t, the number of vertices is v(t) = n*t.
    # The number of possible edges is C(v(t), 2), which is approximately v(t)^2 / 2 for large v(t).
    # So, the rate of edge formation is de/dt = (1/n) * (v(t)^2 / 2).
    # Substituting v(t) = n*t gives de/dt = (1/n) * ((n*t)^2 / 2) = (n * t^2) / 2.
    # Integrating this rate from 0 to c gives the total number of edges E.
    # E = integral from 0 to c of (n * t^2 / 2) dt = n * [t^3 / 6] from 0 to c = n * c^3 / 6.
    print("Step 2: Determine the number of edges (E) at the critical time 'c'.")
    print("The rate of edge formation is de/dt ≈ (1/n) * C(v(t), 2) ≈ (n * t^2) / 2.")
    print("Integrating this rate from time 0 to 'c' gives the total number of edges:")
    print("E = n * c^3 / 6")
    print("-" * 50)

    # Step 3: Apply the condition for the giant component.
    # The giant component emerges when the average degree <k> of the graph equals 1.
    # The average degree is given by <k> = 2 * E / V.
    print("Step 3: Apply the condition for the emergence of the giant component.")
    print("This occurs when the average degree <k> of the graph equals 1.")
    print("<k> = 2 * E / V")
    print("\nSubstituting the expressions for E and V at time 'c':")
    print("<k> = 2 * (n * c^3 / 6) / (n * c)")
    print("Simplifying this expression gives the average degree at time 'c':")
    print("<k> = c^2 / 3")
    print("-" * 50)

    # Step 4: Solve for c.
    # Set the average degree to 1 and solve the resulting equation.
    print("Step 4: Solve for 'c' by setting the average degree to 1.")
    print("This gives us the final equation to solve:")

    # Define and print each number in the final equation: c^2 / 3 = 1
    variable_name = "c"
    power_value = 2
    divisor_value = 3
    result_value = 1
    
    print(f"\n{variable_name}^{power_value} / {divisor_value} = {result_value}\n")

    # Solve the equation
    c_squared = divisor_value * result_value
    c_value = math.sqrt(c_squared)
    
    print("Solving this equation for 'c':")
    print(f"c^2 = {c_squared}")
    print(f"c = sqrt({c_squared})")
    print(f"\nThe exact value of c is the square root of 3, which is approximately {c_value:.4f}.")

if __name__ == '__main__':
    solve_for_c()
