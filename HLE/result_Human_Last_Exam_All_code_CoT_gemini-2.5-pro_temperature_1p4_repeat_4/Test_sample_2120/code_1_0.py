import numpy as np

def solve_average_coordinate():
    """
    Calculates the average value of the complex coordinates z where the
    inverse of the gradient of S(z) approaches zero.
    
    This is equivalent to finding the average value of the coordinates of the poles
    of S'(z), which are the poles of the electric field E(z) and magnetic field B(z).
    """
    
    # Part 1: Poles from the B(z) field, originating from R(z)
    # The denominator of R(z) is the polynomial P1(z) = 4z^4 - z^3 + z^2 + 1.
    # The coefficients are [a_4, a_3, a_2, a_1, a_0]
    p1_coeffs = [4, -1, 1, 0, 1]
    num_roots_p1 = 4
    # According to Vieta's formulas, the sum of roots is -a_{n-1}/a_n
    sum_roots_p1 = -p1_coeffs[1] / p1_coeffs[0]
    
    print(f"Contribution from the first set of B(z) poles:")
    print(f"  Polynomial: {p1_coeffs[0]}z^4 + {p1_coeffs[1]}z^3 + {p1_coeffs[2]}z^2 + {p1_coeffs[3]}z + {p1_coeffs[4]}")
    print(f"  Number of poles: {num_roots_p1}")
    print(f"  Sum of poles: -({p1_coeffs[1]}) / {p1_coeffs[0]} = {sum_roots_p1}")
    print("-" * 20)
    
    # Part 2: Poles from the B(z) field, originating from R(1/z)
    # The denominator of R(1/z) is the polynomial P2(z) = z^4 + z^2 - z + 4.
    p2_coeffs = [1, 0, 1, -1, 4]
    num_roots_p2 = 4
    sum_roots_p2 = -p2_coeffs[1] / p2_coeffs[0]
    
    print(f"Contribution from the second set of B(z) poles:")
    print(f"  Polynomial: {p2_coeffs[0]}z^4 + {p2_coeffs[1]}z^3 + {p2_coeffs[2]}z^2 + {p2_coeffs[3]}z + {p2_coeffs[4]}")
    print(f"  Number of poles: {num_roots_p2}")
    print(f"  Sum of poles: -({p2_coeffs[1]}) / {p2_coeffs[0]} = {sum_roots_p2}")
    print("-" * 20)

    # Part 3: Poles from the E(z) field
    # These are located at z=1 and z=2.
    poles_e = [1, 2]
    num_poles_e = len(poles_e)
    sum_poles_e = sum(poles_e)

    print(f"Contribution from the E(z) poles:")
    print(f"  Pole locations: {poles_e}")
    print(f"  Number of poles: {num_poles_e}")
    print(f"  Sum of poles: {poles_e[0]} + {poles_e[1]} = {sum_poles_e}")
    print("-" * 20)

    # Final Calculation: Average value
    total_num_poles = num_roots_p1 + num_roots_p2 + num_poles_e
    total_sum_poles = sum_roots_p1 + sum_roots_p2 + sum_poles_e
    
    average_value = total_sum_poles / total_num_poles
    
    print("Final Calculation:")
    print(f"  Total number of poles = {num_roots_p1} + {num_roots_p2} + {num_poles_e} = {total_num_poles}")
    print(f"  Total sum of poles = {sum_roots_p1} + {sum_roots_p2} + {sum_poles_e} = {total_sum_poles}")
    print(f"  Average value = Total Sum / Total Number = {total_sum_poles} / {total_num_poles} = {average_value}")

    return average_value

# Execute the function to get the answer
final_answer = solve_average_coordinate()
print(f"\nThe calculated average value is {final_answer}.")
print(f"<<<{final_answer}>>>")
