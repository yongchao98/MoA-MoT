import math

def calculate_braid_index_bound():
    """
    Calculates an upper bound for the braid index of the three-twist knot (6_1)
    using the a-span of its HOMFLYPT polynomial.
    """

    # The HOMFLYPT polynomial for the three-twist knot (6_1) is:
    # P(a, z) = z^4 * a^0  +  a^(-2)*z^4  -  a^(-2)*z^2  -  a^(-4)*z^2
    # We only need the powers of the variable 'a'.
    a_powers = [0, -2, -2, -4]

    # Find the maximum and minimum powers of 'a'
    max_power = max(a_powers)
    min_power = min(a_powers)

    # Calculate the a-span of the polynomial
    a_span = max_power - min_power

    # Apply the Morton-Franks-Williams inequality to find the upper bound
    # b(L) <= span_a / 2 + 1
    # For a knot, the number of components c(L) is 1.
    braid_index_bound = (a_span / 2) + 1

    print("Step 1: The HOMFLYPT polynomial for the three-twist knot (6_1) contains the following powers of 'a':")
    print(sorted(list(set(a_powers))))
    print(f"The maximum power of 'a' is {max_power}.")
    print(f"The minimum power of 'a' is {min_power}.")
    print("")

    print("Step 2: Calculate the a-span of the polynomial.")
    print(f"span_a = max_power - min_power = {max_power} - ({min_power}) = {a_span}")
    print("")

    print("Step 3: Use the Morton-Franks-Williams inequality to find the upper bound on the braid index b(L).")
    print("The formula is: b(L) <= span_a / 2 + 1")
    print("Plugging in the numbers for our final equation:")
    print(f"b(6_1) <= {a_span} / 2 + 1 = {int(a_span / 2)} + 1 = {int(braid_index_bound)}")
    print("")
    
    print(f"An upper bound for the braid index of the three-twist knot is {int(braid_index_bound)}.")

calculate_braid_index_bound()