def calculate_fractional_dehn_twist_coefficient():
    """
    Calculates the fractional Dehn twist coefficient for (D_a D_b)^9
    on a punctured torus.
    """
    # Step 1: Define surface parameters
    g = 1  # genus of the torus
    p = 1  # number of boundary components

    # Step 2: Determine the sequence of twists.
    # The mapping class is (D_a o D_b)^9.
    # The composition means D_b is applied first, then D_a.
    # This sequence is repeated 9 times.
    num_repeats = 9
    base_sequence = ['b', 'a']
    curve_sequence = base_sequence * num_repeats
    k = len(curve_sequence) # Total number of twists, which is 18.

    print(f"The sequence of {k} Dehn twists is along the curves: {curve_sequence}")

    # Step 3: Define the intersection number function
    def intersection_number(c1, c2):
        if c1 == 'a' and c2 == 'b':
            return 1
        elif c1 == 'b' and c2 == 'a':
            return -1
        else:
            return 0

    # Step 4: Calculate the sum of intersection numbers
    # The sum is over all pairs (c_i, c_j) where j < i.
    total_sum = 0
    # i corresponds to the mathematical index k in the formula
    for i in range(k):
        # j corresponds to the mathematical index j in the formula
        for j in range(i):
            ci = curve_sequence[i]
            cj = curve_sequence[j]
            term = intersection_number(ci, cj)
            total_sum += term
            # This part can be uncommented to see every term calculation.
            # if term != 0:
            #    print(f"i={i+1}, j={j+1}: i(c_{i+1}, c_{j+1}) = i({ci}, {cj}) = {term}")


    # Step 5: Calculate the final coefficient
    prefactor = 2 - 2*g - p
    coefficient = total_sum / prefactor

    print("\n--- Calculation Summary ---")
    print("The final equation for the coefficient is: (Sum of intersection numbers) / (2 - 2g - p)")
    print("The numbers in the final equation are:")
    print(f"Sum of intersection numbers: {total_sum}")
    print(f"Denominator (2 - 2*g - p): {prefactor}")
    print(f"Fractional Dehn Twist Coefficient: {total_sum} / {prefactor} = {coefficient}")

calculate_fractional_dehn_twist_coefficient()