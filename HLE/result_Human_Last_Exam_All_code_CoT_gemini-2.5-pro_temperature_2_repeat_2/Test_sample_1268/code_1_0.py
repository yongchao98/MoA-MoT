import math

def display_upper_bound_relation():
    """
    This function displays the upper bound relationship for the maximum norm
    in relation to the covolume for orders in real quadratic fields.

    The relationship is derived from Minkowski's convex body theorem.
    Let k_{k,∞} be the L-infinity norm of the shortest non-zero lattice vector,
    and V be the covolume of the lattice. For a real quadratic field, the dimension n=2.

    Minkowski's theorem, applied with a hypercube (the ball for the L-infinity norm), gives:
    (2 * k_{k,∞})^n <= 2^n * V
    which simplifies to:
    (k_{k,∞})^n <= V

    For n=2, this becomes:
    (k_{k,∞})^2 <= V
    This implies the upper bound for k_{k,∞} is sqrt(V).
    """

    # The key number in the equation is the exponent.
    power = 2

    # Print the derived inequality.
    print("The upper bound for the maximum norm (k_{k,∞}) in relation to the covolume (V) is given by the inequality:")
    print(f"k_{{k,∞}} ^ {power} <= V")
    print("\nThis means that the maximum norm is bounded above by the square root of the covolume:")
    print(f"k_{{k,∞}} <= sqrt(V)")


if __name__ == "__main__":
    display_upper_bound_relation()