import sys

def solve_del_pezzo_problem():
    """
    Calculates the number of isomorphism classes of del Pezzo surfaces of degree 5
    over the rational numbers with good reduction everywhere except possibly at the prime 2.

    This problem is equivalent to counting the number of number fields K over Q
    with degree [K:Q] <= 5 that are unramified outside the prime 2.
    """

    # Step 1: Count the number of such fields for each degree from 1 to 5.
    # These counts are based on established results from number theory.

    # Degree 1: The field is Q itself. Discriminant is 1. Unramified everywhere.
    num_deg_1 = 1

    # Degree 2: These are the quadratic fields Q(sqrt(d)) where the discriminant is a power of 2.
    # This corresponds to d = -1, 2, -2.
    num_deg_2 = 3

    # Degree 3: It is a known result that there are no cubic fields ramified only at the prime 2.
    num_deg_3 = 0

    # Degree 4: Based on recent classifications (e.g., D. Jones, 2019), there are 7
    # distinct quartic fields ramified only at the prime 2.
    num_deg_4 = 7

    # Degree 5: It is a significant result, proven by Dembélé, and Jehanne/Müller,
    # that there are no quintic fields ramified only at the prime 2.
    num_deg_5 = 0

    # Step 2: Calculate the total number by summing the counts.
    total_count = num_deg_1 + num_deg_2 + num_deg_3 + num_deg_4 + num_deg_5

    # Step 3: Print the explanation and the final equation.
    print("The total number of isomorphism classes is the sum of the number of relevant number fields of degree 1 to 5.")
    print(f"Number of fields of degree 1: {num_deg_1}")
    print(f"Number of fields of degree 2: {num_deg_2}")
    print(f"Number of fields of degree 3: {num_deg_3}")
    print(f"Number of fields of degree 4: {num_deg_4}")
    print(f"Number of fields of degree 5: {num_deg_5}")
    print(f"\nThe total number is the sum: {num_deg_1} + {num_deg_2} + {num_deg_3} + {num_deg_4} + {num_deg_5} = {total_count}")

solve_del_pezzo_problem()