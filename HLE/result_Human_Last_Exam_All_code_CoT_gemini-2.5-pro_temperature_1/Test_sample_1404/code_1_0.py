import math

def main():
    """
    Calculates and compares the number of arrangements for two types of candle packages.
    """

    # --- Step 1: Calculate arrangements for the large package ---
    n_large = 9
    # The number of permutations for 9 distinct red candles in 9 positions is 9!
    perms_red_large = math.factorial(n_large)
    # The number of permutations for 9 distinct green candles in 9 positions is 9!
    perms_green_large = math.factorial(n_large)
    # Total arrangements is the product since the placements are independent
    total_arr_large = perms_red_large * perms_green_large

    print("--- Calculation for the Large Package ---")
    print("The large package has 9 distinct red candles for 9 horizontal positions and 9 distinct green for 9 vertical positions.")
    print(f"Number of arrangements for red candles = {n_large}! = {perms_red_large:,}")
    print(f"Number of arrangements for green candles = {n_large}! = {perms_green_large:,}")
    print(f"Total arrangements (Large) = {n_large}! * {n_large}! = {total_arr_large:,}\n")


    # --- Step 2: Calculate arrangements for the small package ---
    n_small = 16
    # The 16 candles (8 red, 8 green, all of unique lengths) are all distinct.
    # The number of permutations for 16 distinct items in 16 positions is 16!
    total_arr_small = math.factorial(n_small)

    print("--- Calculation for the Small Package ---")
    print("The small package has 16 distinct candles in total to be placed in 16 horizontal positions.")
    print(f"Total arrangements (Small) = {n_small}! = {total_arr_small:,}\n")


    # --- Step 3: Compare the number of arrangements ---
    target_ratio = 1260
    # The ratio of small package arrangements to large package arrangements
    calculated_ratio = total_arr_small / total_arr_large

    print("--- Comparison of Arrangements ---")
    print("Is the number of arrangements for the small package 1260 times greater than for the large package?")
    print("We need to calculate the ratio: (Arrangements_Small) / (Arrangements_Large)")
    print("\nThe final equation is:")
    print(f"Ratio = {n_small}! / ({n_large}! * {n_large}!)")
    print(f"Ratio = {total_arr_small:,} / {total_arr_large:,}")
    print(f"Calculated Ratio = {calculated_ratio:.4f}\n")


    # --- Step 4: Conclusion ---
    print("--- Conclusion ---")
    print(f"The calculated ratio is approximately {calculated_ratio:.2f}, which is not equal to the target ratio of {target_ratio}.")
    print("Therefore, the statement that the number of possible arrangements for the small packages is 1260 times greater than for the large packages is false.")


if __name__ == "__main__":
    main()