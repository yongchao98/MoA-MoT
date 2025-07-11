import numpy as np
import matplotlib.pyplot as plt

def get_shape(y1, y2):
    """
    Generates points for the set S given two vectors y1 and y2.
    """
    points = []
    # Generate unit vectors s = (cos(phi), sin(phi))
    for phi in np.linspace(0, 2 * np.pi, 500):
        s = np.array([np.cos(phi), np.sin(phi)])
        # Calculate x = (|dot(y1,s)|^2, |dot(y2,s)|^2)
        x1 = np.dot(y1, s)**2
        x2 = np.dot(y2, s)**2
        points.append([x1, x2])
    return np.array(points)

def plot_shape(ax, points, title):
    """
    Plots the generated points.
    """
    ax.plot(points[:, 0], points[:, 1], 'b-')
    ax.set_title(title)
    ax.set_xlabel(r'$x_1 = |\langle y_1, s \rangle|^2$')
    ax.set_ylabel(r'$x_2 = |\langle y_2, s \rangle|^2$')
    ax.grid(True)
    ax.axis('equal')

# --- Case 1: Orthogonal vectors (should be a simplex) ---
y1_ortho = np.array([1.5, 0])
y2_ortho = np.array([0, 1])
points_ortho = get_shape(y1_ortho, y2_ortho)
# The equation for this simplex is x1/|y1|^2 + x2/|y2|^2 = 1
# x1/2.25 + x2/1 = 1

# --- Case 2: Non-orthogonal vectors (should be an ellipsoid) ---
# y1 = (1, 0), y2 = (cos(pi/3), sin(pi/3))
theta = np.pi / 3
y1_non_ortho = np.array([1.0, 0.0])
y2_non_ortho = np.array([np.cos(theta), np.sin(theta)])
points_non_ortho = get_shape(y1_non_ortho, y2_non_ortho)

# The implicit equation is (x1 + x2 - sin^2(theta))^2 = 4*x1*x2*cos^2(theta)
# Expanding this gives: x1^2 + x2^2 + 2*x1*x2 - 2*(x1+x2)*sin^2(theta) + sin^4(theta) = 4*x1*x2*cos^2(theta)
# Rearranging: x1^2 + x2^2 + 2*x1*x2*(1-2*cos^2(theta)) - 2*(x1+x2)*sin^2(theta) + sin^4(theta) = 0
# A*x1^2 + B*x2^2 + C*x1*x2 + D*x1 + E*x2 + F = 0
A = 1
B = 1
C = 2 * (1 - 2 * np.cos(theta)**2) # 2*cos(2*theta) -> This should be 1-2*cos^2 = -cos(2*theta) oops, let's fix.
# from thought process: x1^2 + x2^2 + 2x1x2(1-2cos^2(theta)) .. -> x1^2+x2^2 - 2x1x2 cos(2theta)... my thought process had a mistake
# Let's re-expand (x1+x2-sin^2(theta))^2 = 4x1x2cos^2(theta)
# x1^2+x2^2+sin^4(t) + 2x1x2 - 2x1sin^2(t) - 2x2sin^2(t) = 4x1x2cos^2(t)
# x1^2+x2^2 + 2x1x2(1-2cos^2(t)) - 2x1sin^2(t) - 2x2sin^2(t) + sin^4(t) = 0
# The cross term is C = 2 * (1 - 2*np.cos(theta)**2) = -2*np.cos(2*theta)
# Re-expanding again: x1^2+x2^2+sin^4 - 2x1sin^2 - 2x2sin^2 + 2x1x2 = 4x1x2cos^2
# x1^2+x2^2+x1x2(2-4cos^2) - 2x1sin^2 - 2x2sin^2 + sin^4 = 0
# C = 2 - 4*np.cos(theta)**2 = 2 - 4 * (1/4) = 1
C_final = 1
D_final = -2*np.sin(theta)**2
E_final = -2*np.sin(theta)**2
F_final = np.sin(theta)**4

print("For the non-orthogonal case with y1=(1,0), y2=(cos(pi/3), sin(pi/3)):")
print("The shape S is an ellipse described by the quadratic equation:")
print("A*x1^2 + B*x2^2 + C*x1*x2 + D*x1 + E*x2 + F = 0")
print(f"where the coefficients are:")
print(f"A = {A}")
print(f"B = {B}")
print(f"C = {C_final}")
print(f"D = {D_final:.4f}")
print(f"E = {E_final:.4f}")
print(f"F = {F_final:.4f}")

# --- Plotting ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
plot_shape(ax1, points_ortho, 'Orthogonal Case (Simplex)')
plot_shape(ax2, points_non_ortho, 'Non-Orthogonal Case (Ellipsoid)')
plt.show()
