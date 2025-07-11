import math

def solve_minimum_sum():
    """
    This function calculates the minimum value of sum(|S_i|) for a given n,
    based on the derived formula.
    
    The problem is to find a general formula for the minimum sum in terms of n.
    The final derived formula is: Minimum Sum = floor(n^2 / 4) + 2.
    
    This script will calculate and display the result for a sample integer n,
    and output the components of the final equation as requested.
    """
    
    # We choose a sample value for n to demonstrate the calculation. Let's use n=12.
    n = 12

    # Calculate the components of the formula.
    n_squared = n**2
    n_squared_by_4 = n_squared / 4
    floor_val = math.floor(n_squared_by_4)
    min_sum = floor_val + 2
    
    print(f"For n = {n}, the minimum value of the sum is calculated using the formula:")
    print(f"Minimum Sum = floor(n^2 / 4) + 2")
    print(f"\nStep-by-step calculation:")
    # Here we output each number in the final equation as requested.
    print(f"1. Square n: {n}^2 = {n_squared}")
    print(f"2. Divide by 4: {n_squared} / 4 = {n_squared_by_4}")
    print(f"3. Take the floor: floor({n_squared_by_4}) = {floor_val}")
    print(f"4. Add 2: {floor_val} + 2 = {min_sum}")
    
    print(f"\nFinal equation for n={n}:")
    print(f"{min_sum} = floor({n}^2 / 4) + 2")


solve_minimum_sum()