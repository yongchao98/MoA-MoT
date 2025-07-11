import math

def solve():
    """
    Calculates the value of phi(n) based on the derived formula.
    """
    try:
        # The problem is stated for n >= 5, but the formula is general.
        n_str = input("Enter the value of n (an integer >= 5): ")
        n = int(n_str)
        if n < 5:
            print("Warning: The problem is defined for n >= 5.")

        # The derivation shows that phi(n) = exp(Tr(P)),
        # where Tr(P) = Tr(Y) - Tr(W).
        # We found Tr(Y) = 2n and Tr(W) = 2/n.
        # So, Tr(P) = 2n - 2/n.

        term1 = 2 * n
        term2 = 2 / n
        trace_p = term1 - term2
        
        # phi(n) = exp(trace_p)
        result = math.exp(trace_p)

        # Output the components of the final equation as requested.
        print(f"The calculation is based on the formula: exp(2*n - 2/n)")
        print(f"For n = {n}:")
        print(f"The first term in the exponent is 2*n = {term1}")
        print(f"The second term in the exponent is 2/n = {term2}")
        print(f"The trace is Tr(P) = {term1} - {term2} = {trace_p}")
        print(f"The final result is phi({n}) = exp({trace_p}) = {result}")

    except ValueError:
        print("Invalid input. Please enter an integer.")
    except Exception as e:
        print(f"An error occurred: {e}")

solve()