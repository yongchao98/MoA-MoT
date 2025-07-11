import sympy

def find_hamiltonicity_threshold():
    """
    This function calculates and explains the d-threshold for Hamiltonicity
    for randomly perturbed graphs with a given minimum degree.
    """
    # Define symbolic variables for n and eta
    n = sympy.Symbol('n')
    eta = sympy.Symbol('eta')

    # Given minimum degree d
    d = n/2 - eta

    print("Step 1: Define the properties of the 'worst-case' graph H_n.")
    print("The d-threshold is determined by the graph H_n that is hardest to make Hamiltonian.")
    print("A canonical example is a graph with a large independent set A.")
    print("The size of the largest possible independent set (alpha) in a graph with minimum degree d is alpha = n - d.")
    
    # Calculate the size of the independent set
    alpha = n - d
    print(f"\nStep 2: Calculate the size of the independent set 'alpha' for d = {d}.")
    print(f"alpha = n - ({d}) = {alpha}")

    print("\nStep 3: Determine the requirement on the random graph G(n,p).")
    print("For H_n U G(n,p) to be Hamiltonian, the random subgraph induced on the independent set A must contain a Hamiltonian path.")
    print("This induced subgraph is a random graph G(alpha, p).")

    print("\nStep 4: State the known threshold for Hamiltonicity in G(m,p).")
    print("The threshold probability for a random graph G(m,p) to have a Hamiltonian path is p_thresh = ln(m)/m.")
    
    # Substitute m with alpha to find the threshold
    p_threshold = sympy.ln(alpha) / alpha
    
    print("\nStep 5: Formulate the final equation for the d-threshold p.")
    print("By substituting m = alpha, we derive the threshold:")

    # Represent the number 2 as an integer in the final output as requested by the prompt.
    p_final_expr = sympy.ln(n/sympy.S(2) + eta) / (n/sympy.S(2) + eta)
    
    print("\nThe final equation for the threshold p is:")
    print("p = log(", end="")
    print("n/", end="")
    print(2, end="")
    print(" + eta) / (n/", end="")
    print(2, end="")
    print(" + eta)")

find_hamiltonicity_threshold()