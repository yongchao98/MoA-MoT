import math

def calculate_max_hausdorff_dimension_of_sidon_set():
    """
    This script explains and calculates the maximum Hausdorff dimension of a Sidon set in [0, 1].

    The logic proceeds in two main steps:
    1.  Derive an upper bound for the dimension.
    2.  Show that this upper bound can be achieved by construction.
    """

    print("This script determines the maximum Hausdorff dimension of a Sidon set in [0,1].\n")

    print("### Step 1: Establishing the Upper Bound")
    print("A set E in [0,1] is a Sidon set if for any x, y, z, w in E, the equation x+y = z+w implies that the pair {x,y} is the same as {z,w}.")
    print("Because E is a Sidon set, the map (x,y) -> x+y is essentially one-to-one.")
    print("This leads to the relationship: dim_H(E+E) = 2 * dim_H(E).")
    print("Since E is a subset of [0,1], its sumset E+E = {x+y | x,y in E} must be a subset of [0,2].")
    print("Any set on the real line has a Hausdorff dimension of at most 1, so dim_H(E+E) <= 1.")
    print("Combining these two facts, we get the inequality:")
    print("2 * dim_H(E) <= 1  =>  dim_H(E) <= 1/2")
    print("This shows the dimension is at most 1/2.\n")


    print("### Step 2: Showing the Bound of 1/2 is Achievable")
    print("We can construct a Sidon set with dimension arbitrarily close to 1/2.")
    print("The construction involves building a self-similar fractal set from a dense *integer* Sidon set.")
    print("An integer Sidon set A is a set of integers where all pairwise sums are unique.")
    print("A key result in number theory states that the maximum size (M) of a Sidon set within {1, ..., N} is approximately sqrt(N).")
    print("The dimension 'd' of the fractal set built from this integer set is given by the formula: d = log(M) / log(N).\n")

    print("### Calculation")
    # Choose a large number N to demonstrate the limit.
    N = 1e12

    # The maximum size M of a Sidon set in {1, ..., N} is M ~ N^(1/2).
    M = N**0.5

    # Calculate the dimension 'd'.
    dimension = math.log(M) / math.log(N)

    print("Let's use a large N to calculate the dimension 'd'.")
    print(f"For N = {N:.0e}, the size of the corresponding integer Sidon set is M ≈ sqrt(N).")
    print(f"M ≈ ({N:.0e})^(1/2) = {M:.0e}")
    print("\nThe dimension 'd' is then calculated using the final equation:")
    print("d = log(M) / log(N)")
    print(f"d = log({M:.0e}) / log({N:.0e})")
    print(f"d = {math.log(M):.4f} / {math.log(N):.4f}")
    print(f"d = {dimension}\n")


    print("As N approaches infinity, this value is exactly 1/2.")
    print("Since the dimension is bounded above by 1/2 and can be made equal to 1/2, the maximum is 1/2.")

if __name__ == '__main__':
    calculate_max_hausdorff_dimension_of_sidon_set()