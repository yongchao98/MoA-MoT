def solve_hat_puzzle():
    """
    This function provides the solution to the hat logic puzzle.
    The reasoning is explained in the text and not derived by this code itself,
    as that would require a complex logical inference engine. This code's
    purpose is to present the final answer clearly.
    """
    
    # The number of people who figured out their hat color
    num_yes = 4
    
    # The specific arrangement of hats around the table
    # B = Black, W = White
    hat_distribution = ['B', 'W', 'W', 'B', 'B', 'W', 'B', 'B', 'W']
    
    # Count the colors to verify
    black_hats = hat_distribution.count('B')
    white_hats = hat_distribution.count('W')
    
    print(f"Total Black Hats: {black_hats}")
    print(f"Total White Hats: {white_hats}")
    print("-" * 20)
    
    print(f"How many people replied 'Yes.'?")
    print(num_yes)
    print("")
    
    print(f"How are black and white hats distributed around the table?")
    # Represent the final configuration as a string for easy reading
    distribution_str = " ".join(hat_distribution)
    print(distribution_str)
    
    # Final check of the key people: the four with white hats
    print("\nThe four people with White hats are the ones who reply 'Yes'.")

solve_hat_puzzle()

# The final answer in the requested format.
# First part is the number of people, second is the distribution string.
final_answer_count = 4
final_answer_distribution = "B W W B B W B B W"