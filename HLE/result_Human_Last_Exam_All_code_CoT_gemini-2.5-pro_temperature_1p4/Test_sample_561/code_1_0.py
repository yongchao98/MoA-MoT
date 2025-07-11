import math

def solve_dimension():
    """
    Calculates the Minkowski-Bouligand dimension for the described fractal keyboard.
    
    The dimension D is solved from the equation: D = log(45) / log(14)
    This is derived from the affinity dimension formula for self-affine fractals:
    N * s1 * s2^(D-1) = 1
    where N=5 (number of black keys), s1=9/14 (y-scaling), and s2=1/14 (x-scaling).
    """
    
    numerator_val = 45
    denominator_val = 14
    
    dimension = math.log(numerator_val) / math.log(denominator_val)
    
    print("The final equation for the dimension is D = log(numerator) / log(denominator)")
    print(f"The number in the numerator of the log is: {numerator_val}")
    print(f"The number in the denominator of the log is: {denominator_val}")
    print(f"The Minkowski-Bouligand dimension of the piano's black keys is: {dimension}")

solve_dimension()