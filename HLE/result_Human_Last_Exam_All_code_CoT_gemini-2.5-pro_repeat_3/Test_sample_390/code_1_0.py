import numpy as np
import matplotlib.pyplot as plt

def solve_and_visualize():
    """
    Analyzes the shape of the set S for a specific 2D case and visualizes it.
    """
    # 1. Define y1 and y2 (already done in the derivation)
    # y1 = np.array([1, 0])
    # y2 = np.array([1/np.sqrt(2), 1/np.sqrt(2)])

    # 2. Generate points for the parameter phi
    phi = np.linspace(0, 2 * np.pi, 400)

    # 3. Calculate x1 and x2 coordinates
    x1 = np.cos(phi)**2
    x2 = 0.5 * (1 + np.sin(2 * phi))

    # 4. Print the derived equation
    center_x, center_y, radius = 0.5, 0.5, 0.5
    print("The shape of the set S in this example is a curve described by the equation of a circle (a special type of ellipse):")
    print(f"(x1 - {center_x})^2 + (x2 - {center_y})^2 = {radius**2}")
    
    # 5. Visualize the result
    plt.figure(figsize=(6, 6))
    plt.plot(x1, x2, label='Set S')
    
    # Plot the full circle for comparison
    circle_phi = np.linspace(0, 2 * np.pi, 100)
    circle_x = center_x + radius * np.cos(circle_phi)
    circle_y = center_y + radius * np.sin(circle_phi)
    plt.plot(circle_x, circle_y, 'r--', label='Equation: $(x_1-0.5)^2 + (x_2-0.5)^2 = 0.25$')
    
    plt.title('Shape of S for a 2D example')
    plt.xlabel('$x_1 = |\langle y_1, s \\rangle|^2$')
    plt.ylabel('$x_2 = |\langle y_2, s \\rangle|^2$')
    plt.grid(True)
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.axis('equal')
    plt.legend()
    # To prevent the program from waiting for the user to close the plot window,
    # we save the figure instead of showing it directly in a blocking way.
    plt.savefig("shape_of_S.png")
    print("\nA plot 'shape_of_S.png' has been generated to visualize the shape.")

solve_and_visualize()
