import sys

def solve_fair_score():
    """
    Calculates the theoretical maximum FAIR compliance score R for a 
    federated knowledge graph system based on given parameters.
    """
    # Problem parameters
    # c: Consistency level of the decentralized identifier resolution
    c = 0.95
    # b: Branching factor of semantic version control
    b = 3

    print(f"Given parameters:")
    print(f"Consistency level (c) = {c}")
    print(f"Branching factor (b) = {b}\n")
    
    # Step 1: Model the individual FAIR metrics based on the parameters
    # Findability (f) is limited by the consistency of identifier resolution.
    f = c
    # Accessibility (a) is also limited by the consistency of identifier resolution.
    a = c
    # Interoperability (i) is inversely affected by the semantic branching factor.
    i = 1 / b
    # Reusability (r) is dependent on both consistent access and semantic clarity.
    r = c * (1 / b)

    print("Derived FAIR metrics:")
    print(f"Findability (f) = c = {f:.4f}")
    print(f"Accessibility (a) = c = {a:.4f}")
    print(f"Interoperability (i) = 1/b = {i:.4f}")
    print(f"Reusability (r) = c * (1/b) = {r:.4f}\n")

    # Step 2: Model the overall FAIR score R as the arithmetic mean of the components
    R = (f + a + i + r) / 4

    # Step 3: Print the final calculation and result
    print("Calculating the overall FAIR score (R) as the average of the metrics:")
    # The prompt requests to output each number in the final equation.
    # We use formatted strings to present the equation clearly.
    equation = f"R = (f + a + i + r) / 4"
    values_equation = f"R = ({f:.4f} + {a:.4f} + {i:.4f} + {r:.4f}) / 4"
    sum_value = f + a + i + r
    result_equation = f"R = {sum_value:.4f} / 4 = {R:.4f}"
    
    print(equation)
    print(values_equation)
    print(result_equation)

    # Final answer as required by the format specification
    # Using sys.stdout to avoid printing a newline character before the answer
    sys.stdout.write(f"\n<<<R = {R:.4f}>>>\n")

solve_fair_score()