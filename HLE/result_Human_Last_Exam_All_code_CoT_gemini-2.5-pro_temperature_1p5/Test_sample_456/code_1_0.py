import numpy as np

def demonstrate_channel_swizzle():
    """
    This function demonstrates the effect of an RGB channel swizzle augmentation.
    It creates a simple 2x2 color image, prints its original state,
    applies a channel swizzle, and prints the result.
    """

    # Create a sample 2x2 image with distinct colors in each quadrant
    # Image has 2 rows, 2 columns, and 3 channels (R, G, B)
    # R: Red, G: Green, B: Blue, Y: Yellow (Red+Green)
    original_image = np.array([
        [[255, 0, 0], [0, 255, 0]],  # Top row: Red, Green
        [[0, 0, 255], [255, 255, 0]]   # Bottom row: Blue, Yellow
    ], dtype=np.uint8)

    print("--- Original Image Pixels (R, G, B) ---")
    print(original_image)
    print("\nThis represents an image with Red, Green, Blue, and Yellow pixels.")

    # --- Perform Channel Swizzle ---
    # The channels are at axis=2. We will shuffle the order of these channels.
    # Original order is [0, 1, 2] -> [R, G, B]
    # Let's define a new order, for example, [2, 0, 1] -> [B, R, G]
    
    # Get the separate channels
    r_channel = original_image[:, :, 0]
    g_channel = original_image[:, :, 1]
    b_channel = original_image[:, :, 2]
    
    # Create the augmented image by stacking channels in a new order: B, R, G
    # We use np.stack to combine them along the third axis.
    swizzled_image = np.stack([b_channel, r_channel, g_channel], axis=2)

    print("\n--- Applying Channel Swizzle (New Order: B, G, R) ---")
    # For demonstration, let's just use permutation. In practice, the permutation would be random.
    # np.random.permutation([0,1,2]) would give a random new order.
    # We will permute from (R,G,B) to (B,R,G)
    
    # Original Red pixel [255, 0, 0] becomes [0, 255, 0] (Green)
    # Original Green pixel [0, 255, 0] becomes [0, 0, 255] (Blue)
    # Original Blue pixel [0, 0, 255] becomes [255, 0, 0] (Red)
    # Original Yellow pixel [255, 255, 0] becomes [0, 255, 255] (Cyan)
    
    permutation_indices = [2, 0, 1] # Corresponds to swapping R->G, G->B, B->R
    augmented_image = original_image[:, :, permutation_indices]

    print("\n--- Augmented Image Pixels after Swizzle ---")
    print(augmented_image)
    print("\nAs you can see, the colors have been re-mapped. This forces the model")
    print("to learn that hue is not a reliable feature, which is critical when")
    print("deploying on grayscale infrared images that have no color information.")


if __name__ == "__main__":
    demonstrate_channel_swizzle()
