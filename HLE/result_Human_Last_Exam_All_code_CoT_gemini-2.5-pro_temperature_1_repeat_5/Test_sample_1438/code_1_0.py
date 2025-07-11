def explain_grassmann_integral():
    """
    Explains the rules of Grassmann variable integration (Berezin integral)
    and their connection to the Pauli exclusion principle.
    """
    print("In the path integral formalism for fermionic systems, the Pauli exclusion principle is enforced by using anti-commuting numbers called Grassmann variables (e.g., η).")
    print("The key property of these variables is that their square is zero (η² = 0), which mirrors the fact that a fermionic state cannot be occupied twice.")
    print("\nIntegration over these variables is defined by the Berezin integral rules.")
    print("For a single Grassmann variable η, the measure 'dη' follows two fundamental rules:")

    # Define the rules
    rule1_lhs = "∫ dη"
    rule1_rhs = 0
    rule2_lhs = "∫ η dη"
    rule2_rhs = 1

    print(f"\n1. Linearity and the integral of a constant:")
    print(f"   The integral of the measure alone is defined to be zero.")
    print(f"   Equation: {rule1_lhs} = {rule1_rhs}")

    print(f"\n2. Normalization:")
    print(f"   This rule provides the non-trivial value and defines the normalization of the integral.")
    print(f"   It is the value of this integral that ensures the formalism is consistent.")
    
    # Print the final equation, explicitly showing the number as requested
    print("\nFinal Defining Equation:")
    print(f"   {rule2_lhs} = {rule2_rhs}")

if __name__ == "__main__":
    explain_grassmann_integral()