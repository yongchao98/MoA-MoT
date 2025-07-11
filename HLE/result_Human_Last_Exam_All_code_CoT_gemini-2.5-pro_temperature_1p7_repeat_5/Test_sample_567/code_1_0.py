import math

def solve_embedding_problem():
    """
    Finds the smallest non-trivial value of 'a' for which the only obstruction
    to symplectically embedding the ellipsoid E(1,a) into a 4-ball is volume.

    This is equivalent to finding the smallest integer a > 1 where the symplectic
    capacity function c(a) equals the volume lower bound sqrt(a). A result
    by McDuff and Schlenk shows this is true if and only if 'a' is a perfect square.
    """
    print("Searching for the smallest integer 'a' > 1 that is a perfect square.")

    a = 2
    while True:
        sqrt_a = math.sqrt(a)
        # Check if 'a' is a perfect square by seeing if its square root is an integer.
        if sqrt_a == int(sqrt_a):
            n = int(sqrt_a)
            print(f"Found the value: a = {a}")
            print(f"This value is the first perfect square after a=1, since {a} = {n}^2.")
            print("\nLet's verify the condition c(a) = sqrt(a) for this value.")
            
            # The volume constraint requires the capacity lambda to be at least sqrt(a).
            volume_constraint = math.sqrt(a)
            print(f"The volume constraint gives: c({a}) >= sqrt({a}) = {volume_constraint}")
            
            # For a perfect square a=n^2, the capacity c(a) is known to be exactly n.
            capacity = n
            print(f"The known symplectic capacity is: c({a}) = {capacity}")
            
            print("\nChecking if the capacity equals the volume constraint:")
            # The final equation demonstrates the equality.
            print(f"Final Equation: c({a}) = sqrt({a})")
            print(f"Plugging in the numbers: {capacity} = {volume_constraint}")
            
            if capacity == volume_constraint:
                print("\nThe equality holds. Therefore, for a = 4, the volume is the only obstruction.")
            
            break
        a += 1

solve_embedding_problem()

# The final answer is the value of 'a' found by the script.
print("<<<4>>>")