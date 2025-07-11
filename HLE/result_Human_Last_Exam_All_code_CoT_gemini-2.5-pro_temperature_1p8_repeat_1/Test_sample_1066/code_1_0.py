import math

def print_lower_bound_formula():
    """
    This function derives and prints the lower bound for the expected watermarking score E[S].

    The expected score for a single token t is E_t = sum_i(p_{t,i} * H_{1/p_{t,i}}),
    where H is the harmonic number.
    The Shannon entropy is H(D_t) = -sum_i(p_{t,i} * log(p_{t,i})).

    A known mathematical inequality provides a lower bound for the difference between these two quantities:
    E_t - H(D_t) >= 1 - (pi^2 / 6)

    Summing over n tokens, we get the total expected score E[S] = sum_t(E_t).
    The average entropy alpha = (1/n) * sum_t(H(D_t)).

    So, E[S] >= sum_t(H(D_t) + 1 - (pi^2 / 6))
               = (sum_t H(D_t)) + n * (1 - pi^2 / 6)
               = n * alpha + n * (1 - pi^2 / 6)

    This script will print the final equation.
    """
    
    n_var = "n"
    alpha_var = "alpha"
    constant_1 = 1
    constant_pi_str = "pi"
    constant_2 = 2
    constant_6 = 6
    
    # Using the math library for the value of pi for a more explicit representation
    pi_val = math.pi
    
    # We are asked to output each number in the final equation.
    # The equation is: E[S] >= n*alpha + n*(1 - pi^2 / 6)
    
    print("The final lower bound for the expected score E[S] is given by the following equation:")
    
    # We construct the string to represent the formula clearly, showing the symbols and numbers.
    equation_str = f"E[S] >= {n_var} * {alpha_var} + {n_var} * ({constant_1} - {constant_pi_str}^({constant_2}) / {constant_6})"
    
    # The print statement below represents each part of the formula, satisfying the prompt's request.
    print(f"The term 'n' is: {n_var}")
    print(f"The term 'alpha' is: {alpha_var}")
    print(f"The number '1' is: {constant_1}")
    print(f"The constant 'pi' ({pi_val:.5f}...) is represented as: {constant_pi_str}")
    print(f"The number '2' (for the exponent) is: {constant_2}")
    print(f"The number '6' (for the divisor) is: {constant_6}")
    
    print("\nResulting equation:")
    print(equation_str)

print_lower_bound_formula()