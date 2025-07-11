import math

def solve_hackenbush_representation():
    """
    Calculates the number of pieces and their colors to represent a
    fraction in red-blue-Hackenbush.
    """
    num = 13
    den = 16

    print(f"Finding the red-blue-Hackenbush representation for the number {num}/{den}.")
    print("-" * 30)

    # The number of pieces 'n' is determined by the denominator, den = 2^n
    if den == 0 or (den & (den - 1)) != 0:
        print("Denominator must be a power of 2.")
        return

    n = int(math.log2(den))

    # We need to solve: num = sum(eps_i * 2^(n-i)) for i=1 to n
    # This is equivalent to finding a signed-digit binary representation of num.
    target = num
    epsilons = []
    
    # We find the coefficients from the most significant to the least significant.
    for i in range(n - 1, -1, -1):
        power_of_2 = 2**i
        
        # The maximum possible sum of the remaining terms.
        # e.g., if we are at 8, the rest (4, 2, 1) can sum to 7.
        remaining_max = power_of_2 - 1

        # Determine the coefficient epsilon for the current power of 2.
        # This algorithm finds the unique non-adjacent form (NAF) representation.
        if target > remaining_max:
            epsilon = 1
        elif target < -remaining_max:
            epsilon = -1
        else: # The target is within the range of the remaining terms
            # We must use the coefficient that makes the new target have the
            # opposite sign or be zero, to ensure the next coefficient can be found.
            # In this unique representation, one must be chosen over the other.
            # A simpler way is to check which choice leaves an even remainder.
            if (target - power_of_2) % 2 == 0:
                epsilon = 1
            else:
                epsilon = -1
    
    # The above logic can be simplified for this problem. A greedy algorithm works.
    target = num
    epsilons = []
    for i in range(n - 1, -1, -1):
        power_of_2 = 2**i
        if abs(target - power_of_2) < abs(target + power_of_2):
            epsilon = 1
        else:
            epsilon = -1
        epsilons.append(epsilon)
        target -= epsilon * power_of_2


    # Build the equation string
    equation_parts = []
    for i, epsilon in enumerate(epsilons):
        term_den = 2**(i + 1)
        if epsilon == 1:
            sign = "+"
            color = "Red"
        else:
            sign = "-"
            color = "Blue"
        
        # Don't show the '+' for the first term
        if i == 0 and sign == "+":
            equation_parts.append(f"1/{term_den}")
        else:
            equation_parts.append(f"{sign} 1/{term_den}")

    equation = f"{num}/{den} = {' '.join(equation_parts)}"
    
    print("The representation is found by solving:")
    print(f"{num} = 8*e1 + 4*e2 + 2*e3 + 1*e4")
    print("\nThe resulting equation is:")
    print(equation)
    print("\nThis corresponds to a stalk of pieces with colors (from ground up):")
    color_sequence = ["Red" if e == 1 else "Blue" for e in epsilons]
    print(color_sequence)
    
    num_pieces = len(epsilons)
    print("\nTotal number of pieces needed is the number of terms in the sum.")
    print(f"Number of pieces: {num_pieces}")


solve_hackenbush_representation()
<<<4>>>