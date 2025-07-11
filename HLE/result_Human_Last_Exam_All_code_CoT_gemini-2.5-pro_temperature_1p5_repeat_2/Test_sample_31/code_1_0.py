import math

def calculate_min_stickers_for_pll():
    """
    Calculates the minimum number of stickers required to identify a PLL case
    using principles of information theory.
    """
    
    # Number of possible PLL cases (states)
    num_cases = 21
    
    # Number of possible colors for a non-top-facing sticker
    num_sticker_colors = 4
    
    # Calculate the information needed in bits to distinguish all cases
    required_bits = math.log2(num_cases)
    
    # Calculate the information provided in bits by observing one sticker
    bits_per_sticker = math.log2(num_sticker_colors)
    
    # The minimum number of stickers is the ceiling of the division
    min_stickers = math.ceil(required_bits / bits_per_sticker)
    
    print("To solve this problem, we determine the minimum information required to distinguish all possible states.")
    print(f"1. Total PLL cases (states) to distinguish: {num_cases}")
    print(f"2. Possible colors per side sticker (outcomes): {num_sticker_colors}\n")
    print("Using information theory:")
    print(f"   - Information needed = log2(number of cases) = log2({num_cases}) = {required_bits:.2f} bits")
    print(f"   - Information from one sticker = log2(number of colors) = log2({num_sticker_colors}) = {bits_per_sticker:.2f} bits\n")
    
    print("The minimum number of stickers required is the ceiling of (Information needed / Information from one sticker).")
    print(f"Final Equation: ceil({required_bits:.2f} / {bits_per_sticker:.2f})")
    print(f"Result: {int(min_stickers)}")
    print("\nTherefore, the minimum number of non-top-facing stickers that must be seen to fully identify the PLL case is 3.")

calculate_min_stickers_for_pll()