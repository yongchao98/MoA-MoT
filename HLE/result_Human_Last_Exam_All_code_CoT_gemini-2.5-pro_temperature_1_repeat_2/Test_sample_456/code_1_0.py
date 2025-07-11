import numpy as np

def demonstrate_horizontal_flip():
    """
    This function demonstrates the effect of a horizontal flip augmentation.
    
    The problem asks to identify the most important augmentation for an animal
    classifier deployed with a fixed side-view camera. Animals can pass in
    front of the camera from left-to-right or right-to-left.
    
    A horizontal flip augmentation is the most critical because it teaches the
    model to recognize an animal regardless of its direction of travel, which is
    a guaranteed source of variation in the deployment scenario.
    
    This code creates a simple image of an arrow pointing right and then
    shows how a horizontal flip creates a realistic new sample of an
    arrow pointing left.
    """
    
    # Create a simple 10x10 "image" with a clear direction (an arrow pointing right)
    # 0 = background, 1 = object
    print("Original Image (Arrow pointing right):")
    original_image = np.array([
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 1, 1, 0, 0, 0, 0],
        [0, 1, 1, 1, 1, 1, 1, 1, 0, 0],
        [0, 1, 1, 1, 1, 1, 1, 1, 0, 0],
        [0, 0, 0, 1, 1, 1, 0, 0, 0, 0],
        [0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    ])
    for row in original_image:
        print(' '.join(map(str, row)))
        
    # Apply a horizontal flip. In numpy, this is done by reversing the order of columns.
    flipped_image = np.fliplr(original_image)
    
    print("\nHorizontally Flipped Image (Arrow pointing left):")
    for row in flipped_image:
        print(' '.join(map(str, row)))

    print("\nThis demonstrates how horizontal flipping creates a valid new training sample for the described scenario.")

if __name__ == "__main__":
    demonstrate_horizontal_flip()
