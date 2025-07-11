import spherogram
import sys

# Suppress snappy low-precision warnings for a cleaner output.
# This does not affect the correctness of the volume calculation (zero vs. non-zero).
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

def solve_knot_proportion():
    """
    This function calculates the proportion of 7-crossing knots that are hyperbolic.
    """
    print("Step 1: Identifying all knot types with 7 crossings.")
    # According to knot theory, there are 7 distinct knot types with 7 crossings.
    # In Alexander-Briggs notation, they are 7_1, 7_2, ..., 7_7.
    # We create a list of their names.
    knot_names = [f'7_{i}' for i in range(1, 8)]
    total_knots = len(knot_names)
    print(f"There are a total of {total_knots} knot types with 7 crossings.")
    print("-" * 30)

    print("Step 2: Checking each knot to see if it is hyperbolic.")
    print("A knot is hyperbolic if the volume of its complement is positive.")
    print("-" * 30)
    
    hyperbolic_knots_count = 0
    hyperbolic_knots_list = []
    
    for name in knot_names:
        # Create a knot object from its name.
        knot = spherogram.Knot(name)
        
        # The volume() method returns the hyperbolic volume of the knot's complement.
        # A positive volume means the knot is hyperbolic.
        # Non-hyperbolic knots (like torus knots) have a volume of 0.
        vol = knot.volume()
        
        if vol > 0:
            hyperbolic_knots_count += 1
            hyperbolic_knots_list.append(name)
            print(f"Knot {name}: Hyperbolic (Volume = {vol:.4f})")
        else:
            print(f"Knot {name}: Non-hyperbolic (Torus Knot) (Volume = {vol:.4f})")
    
    print("-" * 30)
    print("Step 3: Calculating the proportion.")
    print(f"The number of hyperbolic knots found is {hyperbolic_knots_count}.")
    print(f"The total number of knots is {total_knots}.")
    
    proportion = hyperbolic_knots_count / total_knots
    
    print("\nFinal Calculation:")
    print(f"Proportion = (Number of Hyperbolic Knots) / (Total Number of Knots)")
    print(f"Proportion = {hyperbolic_knots_count} / {total_knots}")

if __name__ == '__main__':
    solve_knot_proportion()
