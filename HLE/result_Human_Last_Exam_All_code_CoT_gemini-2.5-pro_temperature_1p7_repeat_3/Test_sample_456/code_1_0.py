import numpy as np

def demonstrate_rgb_swizzle():
    """
    This function demonstrates the RGB channel swizzle augmentation.
    It creates a simple sample image and applies a random permutation
    to its color channels.
    """
    # Create a simple 2x2 image with 4 distinct colored pixels.
    # Each pixel is [R, G, B].
    # The values are high in their respective channels to be clear.
    # Pixel 1: Red, Pixel 2: Green, Pixel 3: Blue, Pixel 4: White
    original_image = np.array([
        [[250, 10, 20],  [10, 250, 20]], # Top row
        [[10, 20, 250],  [250, 250, 250]]  # Bottom row
    ], dtype=np.uint8)

    print("--- Original Image Channels (R, G, B) ---")
    print(original_image)
    print("\nThis represents a 2x2 image. For example, the top-left pixel has a high Red value (250).")

    # The original channel order is [0, 1, 2] corresponding to [R, G, B].
    channel_indices = [0, 1, 2]

    # Randomly shuffle the channel indices.
    # This is the core of the "swizzle" augmentation.
    np.random.shuffle(channel_indices)
    
    # Create the new image by reordering the channels of the original image
    # using the shuffled indices.
    swizzled_image = original_image[:, :, channel_indices]
    
    # Create a mapping from old to new for the print statement
    channel_map = {0: 'R', 1: 'G', 2: 'B'}
    new_order_str = f"[{channel_map[channel_indices[0]]}, {channel_map[channel_indices[1]]}, {channel_map[channel_indices[2]]}]"

    print(f"\n--- Swizzled Image Channels (New order: {new_order_str}) ---")
    print(f"The original channel indices [0, 1, 2] were shuffled to {channel_indices}.")
    print("The original image's channels have been reordered based on this.")
    print(swizzled_image)
    print("\nFor example, if the new order is [G, B, R], the original Red values now appear in the third channel (Blue's original position).")


if __name__ == "__main__":
    demonstrate_rgb_swizzle()