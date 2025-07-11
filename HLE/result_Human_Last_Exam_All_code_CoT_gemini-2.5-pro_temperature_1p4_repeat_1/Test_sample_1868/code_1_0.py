import sys

def calculate_fair_score():
    """
    This script calculates the theoretical maximum FAIR compliance score (R)
    for a federated knowledge graph system based on given system parameters.
    """

    # Given system parameters
    # c: Consistency level of the decentralized identifier resolution
    c = 0.95
    # b: Branching factor of semantic version control
    b = 3

    print("Step 1: Model the maximum values for each FAIR component.")
    print("-" * 50)
    
    # Model f and a based on consistency 'c'
    # The maximum findability and accessibility are capped by the reliability of the identifier system.
    f = c
    a = c
    print(f"The consistency level 'c' is {c}.")
    print(f"This limits Findability (f) and Accessibility (a).")
    print(f"Therefore, maximum f = {f:.2f} and maximum a = {a:.2f}\n")

    # Model i and r based on branching factor 'b'
    # The branching factor introduces semantic ambiguity, which inversely affects interoperability and reusability.
    i = 1 / b
    r = 1 / b
    print(f"The branching factor 'b' is {b}.")
    print(f"This impacts Interoperability (i) and Reusability (r).")
    print(f"Therefore, maximum i = 1/{b} = {i:.3f} and maximum r = 1/{b} = {r:.3f}")
    print("-" * 50)
    
    # Calculate the final FAIR score R
    # R is modeled as the scaled average of the four components.
    R = 10 * (f + a + i + r) / 4

    print("Step 2: Calculate the final FAIR score (R) on a 10-point scale.")
    print("-" * 50)
    print("The score R is calculated as: 10 * (f + a + i + r) / 4")
    print("\nFinal Equation:")
    # sys.stdout.write is used to prevent potential floating point representation issues with f-string formatting
    # when an exact fractional representation is desired in the output string. Here, 1/3 is a repeating decimal.
    sys.stdout.write(f"R = 10 * ({f:.2f} + {a:.2f} + {i:.3f} + {r:.3f}) / 4\n")
    print(f"\nCalculated Result: R = {R:.2f}")

    # Output the final answer in the required format
    print(f"\n<<<{R:.2f}>>>")

# Execute the function
calculate_fair_score()