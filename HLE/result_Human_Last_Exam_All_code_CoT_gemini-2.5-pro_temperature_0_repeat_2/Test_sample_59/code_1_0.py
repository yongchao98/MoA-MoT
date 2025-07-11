import scipy.integrate

def main():
    """
    Calculates the probability of an edge in a jointly exchangeable random graph
    using the graphon model.
    """
    # Define an example graphon function W(u, v).
    # A graphon must be a symmetric function mapping [0,1]x[0,1] to [0,1].
    # Let's use W(u, v) = (u + v) / 2 as an example.
    def W(v, u):  # Note: dblquad integrates with respect to the first arg (v) first
        return (u + v) / 2

    # The probability P(y_ij=1) is the double integral of W(u,v) over the unit square.
    # We use scipy.integrate.dblquad for numerical integration.
    # The arguments are: function, x_min, x_max, y_min, y_max
    # Here, x corresponds to u and y corresponds to v.
    prob, error = scipy.integrate.dblquad(W, 0, 1, 0, 1)

    # The final equation includes the calculated probability.
    # The integral of W(u,v) = (u+v)/2 is:
    # integral from 0 to 1 of [u*v + v^2/2]_0^1 du = integral of (u + 1/2) du
    # = [u^2/2 + u/2]_0^1 = 1/2 + 1/2 = 1. Wait, that's wrong.
    # integral from 0 to 1 ( integral from 0 to 1 of (u+v)/2 dv ) du
    # Inner integral (w.r.t v): [(u*v + v^2/2)/2]_0^1 = (u + 1/2)/2 = u/2 + 1/4
    # Outer integral (w.r.t u): integral from 0 to 1 of (u/2 + 1/4) du
    # = [u^2/4 + u/4]_0^1 = 1/4 + 1/4 = 0.5
    
    lower_bound_u = 0
    upper_bound_u = 1
    lower_bound_v = 0
    upper_bound_v = 1
    
    print("The probability of an edge P(y_ij = 1) is calculated by the formula:")
    print(f"P(y_ij = 1) = integral from u={lower_bound_u} to {upper_bound_u} ( integral from v={lower_bound_v} to {upper_bound_v} W(u,v) dv ) du")
    print("\nFor the example graphon W(u,v) = (u+v)/2, the result is:")
    print(f"P(y_ij = 1) = {prob}")

if __name__ == "__main__":
    main()