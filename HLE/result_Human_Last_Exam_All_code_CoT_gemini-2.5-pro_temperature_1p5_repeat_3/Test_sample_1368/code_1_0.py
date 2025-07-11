import matplotlib.pyplot as plt
import numpy as np

def plot_circle(center, radius, color='black', linestyle='-', label=None, alpha=1.0):
    """Helper function to plot a circle."""
    circle = plt.Circle(center, radius, fill=False, edgecolor=color, linestyle=linestyle, label=label, alpha=alpha)
    ax.add_artist(circle)

def plot_line(p1, p2, color='black', linestyle='-'):
    """Helper function to plot a line segment between two points."""
    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color, linestyle=linestyle)

def plot_infinite_line(p1, p2, color='black', linestyle='-'):
    """Helper function to plot a line that extends to the plot limits."""
    # Calculate the slope and intercept
    if p2[0] - p1[0] == 0:  # Vertical line
        x_vals = [p1[0], p1[0]]
        y_vals = [ax.get_ylim()[0], ax.get_ylim()[1]]
    else:
        slope = (p2[1] - p1[1]) / (p2[0] - p1[0])
        intercept = p1[1] - slope * p1[0]
        x_vals = np.array(ax.get_xlim())
        y_vals = slope * x_vals + intercept
    ax.plot(x_vals, y_vals, color=color, linestyle=linestyle, alpha=0.6)


def get_intersection(c1_center, c1_radius, c2_center, c2_radius):
    """Calculates the intersection points of two circles."""
    x1, y1 = c1_center
    r1 = c1_radius
    x2, y2 = c2_center
    r2 = c2_radius
    
    d = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    
    a = (r1**2 - r2**2 + d**2) / (2 * d)
    h = np.sqrt(r1**2 - a**2)
    
    x_mid = x1 + a * (x2 - x1) / d
    y_mid = y1 + a * (y2 - y1) / d
    
    p3x = x_mid + h * (y2 - y1) / d
    p3y = y_mid - h * (x2 - x1) / d
    
    p4x = x_mid - h * (y2 - y1) / d
    p4y = y_mid + h * (x2 - x1) / d
    
    return (p3x, p3y), (p4x, p4y)

# --- Initial Setup ---
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect('equal')
ax.set_xlim(-3, 3)
ax.set_ylim(-3, 3)
ax.set_title("Constructing an Inscribed Square")
ax.grid(True, linestyle='--', alpha=0.5)

# Given points and circle
P0 = (0, 0)
R = 2.0
P1 = (R, 0)
C0_radius = R

# Plot initial state
plot_circle(P0, C0_radius, color='blue', label='C0 (Given Circle)')
ax.plot(P0[0], P0[1], 'bo', label='P0 (Center)')
ax.plot(P1[0], P1[1], 'ro', label='P1 (Given Point)')

print("Construction sequence: LCCL")
print("--------------------------")

# --- Step 1: Line (L) ---
# Draw a line through P0 and P1 to find P2
P2 = (-R, 0)
plot_infinite_line(P1, P0, color='red', linestyle='--')
ax.plot(P2[0], P2[1], 'ro', label='P2')
print("1. (L)ine: Draw a line through center P0 and point P1. It intersects the circle at a new point, P2.")
print(f"   Vertices found: P1={P1}, P2={P2}")

# --- Step 2: Circle (C) ---
# Draw a circle centered at P1 with radius P1P2
C1_center = P1
C1_radius = np.linalg.norm(np.array(P1) - np.array(P2))
plot_circle(C1_center, C1_radius, color='green', linestyle=':', label='C1')
print("2. (C)ircle: Draw a circle (C1) with center P1 and radius P1-P2.")


# --- Step 3: Circle (C) ---
# Draw a circle centered at P2 with radius P2P1
C2_center = P2
C2_radius = np.linalg.norm(np.array(P2) - np.array(P1))
plot_circle(C2_center, C2_radius, color='orange', linestyle=':', label='C2')
P3, P4 = get_intersection(C1_center, C1_radius, C2_center, C2_radius)
ax.plot(P3[0], P3[1], 'go', label='P3')
ax.plot(P4[0], P4[1], 'go', label='P4')
print(f"3. (C)ircle: Draw a circle (C2) with center P2 and radius P2-P1. The intersections are P3 and P4.")

# --- Step 4: Line (L) ---
# Draw a line through P3 and P4
plot_infinite_line(P3, P4, color='purple', linestyle='--')
# The intersections of this line with the original circle are the other two vertices
V3 = (0, R)
V4 = (0, -R)
ax.plot(V3[0], V3[1], 'mo', label='V3')
ax.plot(V4[0], V4[1], 'mo', label='V4')
print("4. (L)ine: Draw a line through P3 and P4. It intersects the original circle at V3 and V4.")
print(f"   Vertices found: V3=({V3[0]:.2f}, {V3[1]:.2f}), V4=({V4[0]:.2f}, {V4[1]:.2f})")


# --- Final Result ---
# Draw the final square
square_vertices = [P1, V3, P2, V4, P1]
square_x, square_y = zip(*square_vertices)
ax.plot(square_x, square_y, color='black', linewidth=2, label='Final Square')

print("\nFinal square vertices:")
final_vertices = [P1, V3, P2, V4]
for i, v in enumerate(final_vertices, 1):
    print(f"Vertex {i}: ({v[0]:.2f}, {v[1]:.2f})")


ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.show()

# Final Answer as requested by the format
final_answer = 'LCCL'