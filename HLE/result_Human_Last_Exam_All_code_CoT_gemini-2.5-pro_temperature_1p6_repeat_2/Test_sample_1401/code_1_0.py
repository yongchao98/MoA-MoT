import math

def get_lower_bound_explanation():
    """
    This function provides the theoretical lower bound for the number of queries
    in the specified SQ learning problem.
    """
    
    # Define the components of the theoretical lower bound expression.
    # The problem asks for the minimum number of queries, denoted by Q.
    # The dimension of the input space is d.
    # The lower bound is a quasi-polynomial function of d.
    
    queries_variable = "Q"
    inequality_symbol = ">="
    base_of_exponent = "exp"
    asymptotic_notation = "Ω" # This is Big-Omega, representing a lower bound.
    exponent_expression = "log(d)^2"

    # Assemble the final equation string.
    # The result is Q >= exp(Ω(log(d)^2)).
    # This means the number of queries must be at least on the order of exp(log(d)^2).
    equation = f"{queries_variable} {inequality_symbol} {base_of_exponent}({asymptotic_notation}({exponent_expression}))"

    print("For any SQ algorithm learning a poly(d)-sized two-hidden-layer ReLU network")
    print("over N(0, I_d) to a squared loss of 1/poly(d) with a non-negligible")
    print("query tolerance, the minimum number of queries (Q) required is given by the")
    print("following theoretical lower bound:")
    print("\n" + equation + "\n")
    print("This can also be written as: Q >= d^(Ω(log(d)))")
    print("This is a quasi-polynomial lower bound, indicating the problem is computationally hard for the SQ model.")


if __name__ == "__main__":
    get_lower_bound_explanation()
