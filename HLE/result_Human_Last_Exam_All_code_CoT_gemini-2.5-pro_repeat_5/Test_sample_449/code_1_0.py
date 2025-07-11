import math

def solve():
    """
    Calculates the probability that the conditioned random walk never enters the
    set of the four neighbors of the origin.
    """
    # Starting point and its norm
    x0 = (3000, 4000)
    norm_x0 = math.sqrt(x0[0]**2 + x0[1]**2)

    # The set of neighbors of the origin
    y = (1, 0) # Representative neighbor

    # The potential kernel at the neighbor y
    a_y = 1.0

    # Constants for the asymptotic formula of the potential kernel
    gamma = 0.5772156649015328
    ln2 = math.log(2)
    pi = math.pi
    
    # Asymptotic formula for a(x0):
    # a(x0) ~ (2/pi) * ln|x0| + (2*gamma + 3*ln2)/pi
    ln_norm_x0 = math.log(norm_x0)
    
    # C = (2*gamma + 3*ln2)/pi
    constant_term = (2 * gamma + 3 * ln2) / pi
    
    # a(x0)
    a_x0 = (2 / pi) * ln_norm_x0 + constant_term

    # The probability of hitting the set N is v(x0) = a_y / a(x0)
    prob_hit = a_y / a_x0

    # The probability of never hitting the set N is 1 - v(x0)
    prob_avoid = 1 - prob_hit

    print("The probability is calculated using the formula: P = 1 - a(y) / a(x0)")
    print(f"where y = {y} is a neighbor of the origin, and x0 = {x0}.")
    print("-" * 30)
    print("Values used in the calculation:")
    print(f"a(y) = a(1,0) = {a_y}")
    print(f"|x0| = {norm_x0}")
    print(f"gamma = {gamma:.6f}")
    print(f"ln(2) = {ln2:.6f}")
    print(f"pi = {pi:.6f}")
    print(f"ln(|x0|) = ln({norm_x0}) = {ln_norm_x0:.6f}")
    print(f"The potential kernel a(x0) is approximated as: a(x0) \u2248 {a_x0:.6f}")
    print("-" * 30)
    print("Final equation:")
    print(f"P = 1 - {a_y} / {a_x0:.6f}")
    print(f"P = 1 - {prob_hit:.6f}")
    print(f"P = {prob_avoid:.6f}")
    print("-" * 30)
    print(f"The approximate answer with two significant digits is: {prob_avoid:.2f}")

solve()