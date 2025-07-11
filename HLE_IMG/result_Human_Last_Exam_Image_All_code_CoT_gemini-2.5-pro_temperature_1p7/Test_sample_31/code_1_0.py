def calculate_rasmussen_invariant_T25():
    """
    Calculates the Rasmussen invariant for the knot in the image.
    The knot is identified as the T(2,5) torus knot.
    The formula for the Rasmussen invariant of a positive torus knot T(p,q)
    is s = (p-1)*(q-1).
    """

    # Parameters for the T(2,5) torus knot
    p = 2
    q = 5

    print("The knot shown in the image is the positive torus knot T(p,q) with p=2 and q=5.")
    print("The formula for the Rasmussen s-invariant for a positive torus knot is:")
    print("s = (p - 1) * (q - 1)")
    print("\nCalculating for the T(2,5) knot:")

    # Perform the calculation step by step
    p_minus_1 = p - 1
    q_minus_1 = q - 1
    invariant = p_minus_1 * q_minus_1

    # Print the equation with the values substituted
    print(f"s = ({p} - 1) * ({q} - 1)")
    print(f"s = {p_minus_1} * {q_minus_1}")
    print(f"s = {invariant}")
    print(f"\nThe Rasmussen invariant of the knot is {invariant}.")

# Execute the function to print the result
calculate_rasmussen_invariant_T25()