import scipy.integrate

# A jointly exchangeable random graph can be represented by a graphon W(u,v).
# A graphon is a symmetric, measurable function from [0,1]x[0,1] to [0,1].
# The probability of an edge y_ij is given by P(y_ij=1) = integral(W(u,v) du dv) over the unit square.

# We define an example graphon W(u, v) = u * v to demonstrate the calculation.
# Note that for dblquad, the order of arguments in the function should be (y, x) -> (v, u).
def graphon_W(v, u):
    """
    Example graphon function W(u,v) = u * v.
    It is symmetric: W(u,v) = W(v,u).
    It maps [0,1]x[0,1] to [0,1].
    """
    return u * v

# We calculate the double integral of the graphon over the unit square [0,1] x [0,1].
# The limits for the inner integral (du) are 0 to 1.
# The limits for the outer integral (dv) are 0 to 1.
# dblquad returns a tuple (result, error_estimate).
prob, err = scipy.integrate.dblquad(graphon_W, 0, 1, lambda u: 0, lambda u: 1)

# The result of the integral is the unconditional probability of an edge.
print("For a jointly exchangeable graph defined by the example graphon W(u,v) = u * v:")

# The final equation includes the definition and the calculated value.
# The numbers in the equation are the limits of integration (0 and 1) and the final probability.
print(f"P(y_ij=1) = Integral from 0 to 1 of (Integral from 0 to 1 of (u * v) du) dv = {prob}")
