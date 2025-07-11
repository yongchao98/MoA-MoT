import math

def calculate_arrangements():
    """
    Calculates the number of arrangements for both large and small packages
    and determines if the small package has 1260 times more arrangements
    than the large one.
    """
    
    # Step 1: Calculate arrangements for the large package
    # There are 9 distinct red candles for 9 horizontal positions (9!) and
    # 9 distinct green candles for 9 vertical positions (9!).
    # Total arrangements = 9! * 9!
    n_large_red = 9
    n_large_green = 9
    
    fact_9 = math.factorial(n_large_red)
    
    num_arrangements_large = fact_9 * fact_9

    print("--- Large Package Calculation ---")
    print(f"The large package has {n_large_red} red candles for {n_large_red} horizontal positions, and {n_large_green} green for {n_large_green} vertical positions.")
    print("The number of arrangements is the product of the permutations for each color.")
    print(f"Arrangements = {n_large_red}! * {n_large_green}!")
    print(f"Result = {fact_9} * {fact_9} = {num_arrangements_large}")
    print("-" * 35)

    # Step 2: Calculate arrangements for the small package
    # There are 8 red and 8 green distinct candles, for a total of 16 distinct candles.
    # All 16 must be placed in 16 horizontal positions.
    # Total arrangements = 16!
    n_small_total = 16
    num_arrangements_small = math.factorial(n_small_total)
    
    print("\n--- Small Package Calculation ---")
    print(f"The small package has a total of {n_small_total} distinct candles to be placed in {n_small_total} positions.")
    print("The number of arrangements is the permutation of all 16 candles.")
    print(f"Arrangements = {n_small_total}!")
    print(f"Result = {num_arrangements_small}")
    print("-" * 35)

    # Step 3: Calculate the ratio and check the condition
    ratio = num_arrangements_small / num_arrangements_large
    
    print("\n--- Comparison ---")
    print("The problem asks if the number of arrangements for the small package is 1260 times greater than for the large package.")
    print(f"Ratio = (Small Package Arrangements) / (Large Package Arrangements)")
    print(f"Ratio = {num_arrangements_small} / {num_arrangements_large} = {ratio:.4f}")
    
    is_true = (ratio == 1260)
    
    print(f"\nIs the calculated ratio equal to 1260? {is_true}")


if __name__ == "__main__":
    calculate_arrangements()
    print("\n<<<No>>>")