import sympy

def solve_whitening_filter():
    """
    This function calculates the coefficients for a simple first-order FIR whitening filter
    W(D) = w0 + w1*D that aims to cancel the first anti-causal term of the output.
    """
    # Define symbols for the filter coefficients
    w0, w1 = sympy.symbols('w0 w1')

    # Given channel coefficients q_k
    # q_k = 2 - (|k|-1)/2 for k odd
    # q_k = 5/3 - |k|/2 for k even
    q_minus_1 = 2 - (abs(-1) - 1) / 2
    q_minus_2 = 5/3 - abs(-2) / 2

    # The condition for the first anti-causal term h_-1 to be zero is:
    # h_-1 = w0*q_-1 + w1*q_-2 = 0
    equation = sympy.Eq(w0 * q_minus_1 + w1 * q_minus_2, 0)

    # We can choose w0 = 1 to find a specific solution for w1
    solution = sympy.solve(equation.subs(w0, 1), w1)
    
    w0_val = 1
    w1_val = solution[0]

    # Print the resulting filter equation
    print("Based on simplifying the problem to finding a first-order FIR filter W(D) = w0 + w1*D")
    print("that cancels the h_{-1} term, we get the following equation:")
    print(f"h_{-1} = {w0_val}*({q_minus_1}) + w1*({q_minus_2}) = 0")
    print(f"Solving for w1 with w0=1 gives w1 = {w1_val}.")
    print("\nThe resulting filter is:")
    # The final equation needs to have each number outputted.
    print(f"W(D) = {w0_val} + ({w1_val})*D")
    print(f"W(D) = {w0_val} - {abs(w1_val)}*D")


solve_whitening_filter()