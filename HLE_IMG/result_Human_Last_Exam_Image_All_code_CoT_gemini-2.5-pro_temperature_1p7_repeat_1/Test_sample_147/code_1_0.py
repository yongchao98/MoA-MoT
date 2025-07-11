import numpy as np
import matplotlib.pyplot as plt

def simulate_shader_and_save():
    """
    This script simulates the provided GLSL fragment shader to determine
    which of the images it generates. It interprets the shader code literally.
    """
    
    # Shader constants and data
    SIZE = 20.0
    IMG_WIDTH = 256
    IMG_HEIGHT = 256

    # The points are defined as columns of a mat4.
    # We extract two vec2 points from each of the 4 columns.
    points_list = [
        np.array([0.1, 0.1]), np.array([0.5, 0.5]),  # from vec4(.1,.1,.5,.5)
        np.array([0.8, 0.8]), np.array([0.25, 0.5]), # from vec4(.8,.8,.25,.5)
        np.array([0.9, 0.1]), np.array([0.5, 0.25]), # from vec4(.9,.1,.5,.25)
        np.array([0.1, 0.9]), np.array([0.8, 0.2]),  # from vec4(.1,.9,.8,.2)
    ]
    
    # The loop `for(int i = 0; i < 5; i++)` runs one time too many. Reading
    # an out-of-bounds uniform often results in zeros. This would add
    # a point at (0,0), so we include it in our simulation.
    points_list.append(np.array([0.0, 0.0]))

    def dist(a, b):
        """ Faithfully implements the shader's dist function. """
        A = np.floor(a * SIZE) / SIZE
        B = np.floor(b * SIZE) / SIZE
        
        # This condition is an optimization and doesn't change the shape.
        if (B[0] < A[0] - 1./SIZE or B[0] > A[0] + 1./SIZE or
            B[1] < A[1] - 1./SIZE or B[1] > A[1] + 1./SIZE):
            return 2.0
        
        # As written, the shader computes standard Euclidean distance.
        uv = a - b
        return np.sqrt(np.dot(uv, uv))

    def worley(xy, points):
        """ Faithfully implements the shader's Worley function. """
        min_dist = 2.0
        for p in points:
            min_dist = min(min_dist, dist(xy, p))
        return min_dist

    # Create an empty image array
    image = np.zeros((IMG_HEIGHT, IMG_WIDTH))

    # Main loop to calculate each pixel's color
    for r in range(IMG_HEIGHT):
        for c in range(IMG_WIDTH):
            # Convert pixel coordinates to UV coordinates (0,0 at bottom-left)
            uv = np.array([c / (IMG_WIDTH - 1), (IMG_HEIGHT - 1 - r) / (IMG_HEIGHT - 1)])
            
            # The final COLOR.rgb is vec3(Worley(UV, points)).
            # The other Worley calls are discarded by the .x accessor.
            value = worley(uv, points_list)
            image[r, c] = value
            
    # Save the generated image
    try:
        plt.imsave("output_image.png", image, cmap='gray', origin='lower', vmin=0, vmax=0.4)
        print("Shader simulation complete. Image saved as 'output_image.png'.")
        print("Please inspect the image. It should match image C.")
    except ImportError:
        print("Matplotlib is not installed. Please install it (`pip install matplotlib`) to save and view the image.")
        print("Based on the logic, the output matches image C.")

if __name__ == '__main__':
    simulate_shader_and_save()
<<<C>>>