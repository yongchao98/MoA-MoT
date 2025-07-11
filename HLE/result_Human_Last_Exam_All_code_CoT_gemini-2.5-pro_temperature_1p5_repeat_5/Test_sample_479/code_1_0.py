import sys

def solve_genus_problem():
    """
    Solves for the maximal genus of a surface based on a geometric theorem.

    The problem states that Ω is a compact region in R^3 with a smooth connected
    boundary ∂Ω, and the mean curvature vector of ∂Ω never vanishes.

    This implies the mean curvature scalar H is either always positive or always
    negative. A fundamental theorem in differential geometry (by Aleksandrov and Hopf)
    states that such a surface must be topologically equivalent to a sphere.

    We can now determine the genus of a sphere using the formula relating
    the Euler characteristic (χ) and the genus (g):
    χ = 2 - 2g

    For a sphere, the Euler characteristic is 2. We can substitute this into
    the formula and solve for g.
    """
    # Known value for the Euler characteristic of a sphere.
    euler_characteristic_chi = 2

    # The formula is chi = 2 - 2*g. We print out the numbers in the final equation.
    constant_term = 2
    genus_coefficient = 2
    
    print("Based on the Aleksandrov-Hopf theorem, the surface must be a sphere.")
    print("The Euler characteristic (χ) of a sphere is 2.")
    print("The formula relating Euler characteristic (χ) and genus (g) is: χ = 2 - 2g")
    print("\nSubstituting χ = 2 into the formula, we get the equation:")
    # Using python's f-string to output the numbers in the equation
    print(f"{euler_characteristic_chi} = {constant_term} - {genus_coefficient}*g")
    
    print("\nSolving the equation for g:")
    step1_result = euler_characteristic_chi - constant_term
    print(f"1. Subtract {constant_term} from both sides: {euler_characteristic_chi} - {constant_term} = -{genus_coefficient}*g")
    print(f"   {step1_result} = -{genus_coefficient}*g")
    
    final_genus_g = step1_result / -genus_coefficient
    print(f"2. Divide by -{genus_coefficient}: g = {step1_result} / -{genus_coefficient}")
    print(f"   g = {int(final_genus_g)}")

    print("\nSince the only possible genus is 0, the maximal genus is 0.")


solve_genus_problem()