import sys

def solve_tori_count():
    """
    This function calculates the number of Fq-rational maximal tori for a reductive group
    of type E_8 and prints a detailed explanation of the solution.
    """

    # The reductive group G is of type E_8.
    group_type = "E_8"

    # The degrees of the fundamental invariants for the Weyl group of type E_8.
    degrees = [2, 8, 12, 14, 18, 20, 24, 30]

    # The number of positive roots N is calculated by the formula: sum(d_i - 1)
    # where d_i are the degrees of the fundamental invariants.
    num_positive_roots = sum(d - 1 for d in degrees)

    # Explanation of the method
    print("To find the exact number of Fq-rational maximal tori of a reductive group G of type E_8,")
    print("we rely on a theorem from the theory of algebraic groups over finite fields.")
    print("\nTheorem (Steinberg): The number of Fq-rational maximal tori in a connected reductive group")
    print("over the finite field Fq is given by q^N, where N is the number of positive roots")
    print("of the associated root system.")
    
    # Calculation steps
    print(f"\nFor a group of type {group_type}, the number of positive roots N can be found using the")
    print("degrees (d_i) of the fundamental invariants of its Weyl group. The formula is:")
    print("N = Sum(d_i - 1)")
    
    print(f"\nThe degrees for {group_type} are: {degrees}")

    print("\nWe calculate N with each number in the equation:")
    
    # Building the string for the equation with each number explicitly shown
    calc_expression = " + ".join([f"({d}-1)" for d in degrees])
    calc_values = " + ".join([str(d - 1) for d in degrees])
    
    print(f"N = {calc_expression}")
    print(f"N = {calc_values}")
    print(f"N = {num_positive_roots}")

    # Final result
    print("\nWith N calculated, the number of Fq-rational maximal tori is q^N.")
    print(f"Therefore, the exact number is q^{num_positive_roots}.")

# Execute the function to print the solution.
solve_tori_count()

# Redirecting stderr to stdout for the final answer format if needed,
# though here print will suffice.
# The final answer in the required format
# We use a trick to make it part of the output without extra print statements from the environment
sys.stdout.write("\n<<<q^120>>>\n")