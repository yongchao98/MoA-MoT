import scipy.integrate as integrate

def calculate_edge_probability():
    """
    Calculates the probability of an edge in a jointly exchangeable random graph
    for a given example graphon.
    """
    # The graphon function W(u,v) defines the probability of an edge given
    # the latent variables u and v for the two nodes.
    # We use the example W(u,v) = u * v.
    def W(u, v):
        return u * v

    # The probability of an edge is the integral of W(u,v) over the unit square [0,1]x[0,1].
    lower_limit = 0
    upper_limit = 1

    # Use scipy's double quadrature function to compute the integral.
    # It returns the result and an estimate of the absolute error.
    prob, error = integrate.dblquad(W, lower_limit, upper_limit, lambda x: lower_limit, lambda x: upper_limit)

    # Print the explanation and the final equation with the calculated values.
    print("The probability P(y_ij = 1) is given by the integral of the graphon W(u,v).")
    print("For the example graphon W(u,v) = u * v, the probability is calculated as follows:")
    print(f"P(y_ij = 1) = integral from {lower_limit} to {upper_limit} ( integral from {lower_limit} to {upper_limit} (u * v) du ) dv = {prob:.4f}")
    print(f"(Numerical integration error estimate: {error:.2e})")

if __name__ == '__main__':
    calculate_edge_probability()