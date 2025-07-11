import math

def calculate_virial_prefactor(n):
    """
    Calculates the system-independent prefactor c_n for the fully f-connected
    diagram in the virial expansion.

    The formula for the prefactor is c_n = -(n-1) / n!

    Args:
        n (int): The number of nodes in the diagram (must be an integer >= 2).
    """
    if not isinstance(n, int) or n < 2:
        print("Error: The number of nodes 'n' must be an integer greater than or equal to 2.")
        return

    # Calculate the numerator and denominator of the fraction
    numerator = n - 1
    denominator = math.factorial(n)
    
    # Calculate the value of c_n
    c_n_value = -numerator / denominator

    # Print the results, showing each number in the equation
    print(f"Calculation for the prefactor c_n for n = {n}:")
    print(f"The general formula is: c_n = -(n-1) / n!")
    print(f"Plugging in n = {n}:")
    print(f"  Numerator (n-1) = {n} - 1 = {numerator}")
    print(f"  Denominator (n!) = {n}! = {denominator}")
    print(f"The final equation is: c_{n} = -{numerator} / {denominator}")
    print(f"The numerical value is: c_{n} = {c_n_value}")

if __name__ == '__main__':
    # You can change this value to calculate the prefactor for a different 'n'
    n_value = 5
    calculate_virial_prefactor(n_value)