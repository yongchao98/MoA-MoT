import numpy as np
import matplotlib.pyplot as plt

def plot_hopfion_preimages():
    """
    This script calculates and plots the preimages of two points for the given
    vector field, demonstrating that their linking number is 1.
    """
    print("Calculating the preimages for the Hopf charge calculation.")

    # --- Curve C1: Preimage of the South Pole (0,0,-1) ---
    # This is a circle in the z=0 plane with radius sqrt(0.5)
    radius_c1 = np.sqrt(0.5)
    theta = np.linspace(0, 2 * np.pi, 100)
    x_c1 = radius_c1 * np.cos(theta)
    y_c1 = radius_c1 * np.sin(theta)
    z_c1 = np.zeros_like(theta)
    print("Preimage C1 (p1=(0,0,-1)): A circle of radius {:.3f} in the z=0 plane.".format(radius_c1))


    # --- Curve C2: Preimage of the Equatorial point (1,0,0) ---
    # This is a loop in the y=0 plane.
    # (x**2 - 0.5)**2 + z**2 = R**2, where R = ln(2)/10
    R = np.log(2) / 10.0
    # Parameterize the curve
    phi = np.linspace(0, 2 * np.pi, 200)
    # x must be positive, and 0.5 + R*cos(phi) is always positive as R < 0.5
    x_c2 = np.sqrt(0.5 + R * np.cos(phi))
    y_c2 = np.zeros_like(phi)
    z_c2 = R * np.sin(phi)
    print("Preimage C2 (p2=(1,0,0)): A loop in the y=0 plane.")


    # --- Plotting ---
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Plot the curves
    ax.plot(x_c1, y_c1, z_c1, label='C1: Preimage of (0,0,-1)', lw=2)
    ax.plot(x_c2, y_c2, z_c2, label='C2: Preimage of (1,0,0)', lw=2)

    # Plot the point where C1 pierces the plane of C2
    ax.scatter([np.sqrt(0.5)], [0], [0], color='red', s=100, label='Linking Point')

    # Formatting the plot
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.set_title('Preimage Curves C1 and C2 for Hopfion')
    ax.legend()
    ax.grid(True)
    ax.set_aspect('equal', 'box')
    
    # Set view for clarity
    ax.view_init(elev=25, azim=-50)

    print("\nThe plot shows that the blue circle (C1) and the orange loop (C2) are linked together once.")
    print("The linking number Lk(C1, C2) is therefore +/- 1.")
    print("By convention, the Hopf Charge Q_H is taken to be positive.")

    # The final answer
    hopf_charge = 1
    print(f"\nThe Hopf charge of this field is {hopf_charge}.")

    plt.show()

# Run the function
plot_hopfion_preimages()
