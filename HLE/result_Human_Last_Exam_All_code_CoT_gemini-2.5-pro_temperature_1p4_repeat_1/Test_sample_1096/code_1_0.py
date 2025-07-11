import math

def solve_optimal_steps():
    """
    This function solves for the optimal two-step gradient descent learning rates
    based on the problem description.

    The problem is interpreted as a riddle where the condition number M is
    encoded in the definition of S.
    The equation M^2 + 6M + 1 = 2M^2 - 2M + 1 yields M=8.
    """
    M = 8.0
    mu = 1.0

    # From the general formula for optimal step sizes:
    # gamma = (4*(M+mu) +/- 2*sqrt(2)*(M-mu)) / (M^2 + 6*M*mu + mu^2)
    
    # Numerator parts
    term1_num = 4 * (M + mu)
    term2_num = 2 * math.sqrt(2) * (M - mu)
    
    # Denominator
    denominator = M**2 + 6*M*mu + mu**2
    
    # The final expression for the step sizes is of the form:
    # (a +/- b*sqrt(c)) / d
    # We are asked to output each number in the final equation.
    a = 4 * (M + mu)
    b = 2 * (M - mu) # we have b * sqrt(2)
    c = 2
    d = M**2 + 6 * M * mu + mu**2
    
    # Simplified forms for printing
    a_int = int(a)
    b_int = int(b/math.sqrt(c))
    c_int = int(c)
    d_int = int(d)

    print("The best choice for the pair (gamma_1, gamma_2) is given by the expression:")
    print(f"(a +/- b*sqrt(c)) / d")
    print("\nThe numbers in this final equation are:")
    print(f"a = {a_int}")
    print(f"b = {b_int}")
    print(f"c = {c_int}")
    print(f"d = {d_int}")
    
    gamma1_val = (term1_num - term2_num) / denominator
    gamma2_val = (term1_num + term2_num) / denominator
    
    print("\nThe two optimal step sizes are therefore:")
    print(f"gamma_1 = ({a_int} - {b_int}*sqrt({c_int})) / {d_int} approx {gamma1_val}")
    print(f"gamma_2 = ({a_int} + {b_int}*sqrt({c_int})) / {d_int} approx {gamma2_val}")

if __name__ == '__main__':
    solve_optimal_steps()