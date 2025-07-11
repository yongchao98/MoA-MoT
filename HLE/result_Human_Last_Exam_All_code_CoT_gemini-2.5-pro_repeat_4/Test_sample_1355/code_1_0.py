import math

def calculate_conductance_moment(n):
    """
    Calculates the n-th statistical moment of the dimensionless conductance 'g'
    for a disordered Majorana wire at the critical point (symmetry class D).
    The formula is: <g^n> = Gamma(n + 0.5) / (sqrt(pi) * Gamma(n + 1))
    
    Args:
        n (int): The order of the moment.
        
    Returns:
        float: The value of the n-th moment.
    """
    numerator = math.gamma(n + 0.5)
    denominator = math.sqrt(math.pi) * math.gamma(n + 1)
    return numerator / denominator

def main():
    """
    Calculates the ratio of the fourth moment of conductance to its average value.
    """
    # Calculate the average value (1st moment, n=1)
    avg_conductance = calculate_conductance_moment(1)

    # Calculate the fourth moment (n=4)
    fourth_moment = calculate_conductance_moment(4)

    # Calculate the ratio
    if avg_conductance == 0:
        ratio = float('inf')
    else:
        ratio = fourth_moment / avg_conductance

    print("This script calculates the ratio of the fourth statistical moment of the dimensionless conductance to its average value for a disordered Majorana wire at the critical point.")
    print("-" * 50)
    print(f"The average conductance <g> is: {avg_conductance}")
    print(f"The fourth moment <g^4> is: {fourth_moment}")
    print("-" * 50)
    print("The final ratio is calculated as follows:")
    print(f"<g^4> / <g> = {fourth_moment} / {avg_conductance} = {ratio}")
    
    # The exact result is 35/64
    print(f"\nThe exact value of the ratio is 35/64, which is {35/64}.")

if __name__ == "__main__":
    main()