import math

def solve():
    """
    Searches for a box with smaller surface area than the original 12x12x12 cube
    that can still hold 27 balls under the given packing constraints.
    """
    original_sa = 6 * (12 * 12)
    best_box = None
    min_sa = original_sa

    # The most compact configuration of 27 balls we can construct requires a
    # bounding box for their centers of 10x10x6 cm. We need to find a container
    # that can fit this cluster.
    # The space available for centers inside a container L x W x H is (L-4) x (W-4) x (H-4).
    # We check if a 10x10x6 object can fit inside an (L-4)x(W-4)x(H-4) space.
    cluster_dims = sorted([10, 10, 6])

    # We only need to check integer dimensions up to a reasonable limit.
    # If a dimension is too large, the surface area will also be large.
    # We can establish a rough upper bound, but searching up to 20 is sufficient.
    limit = 20
    for l in range(1, limit + 1):
        for w in range(1, l + 1):
            for h in range(1, w + 1):
                # Calculate surface area
                sa = 2 * (l * w + l * h + w * h)

                # We are looking for a more efficient box
                if sa < min_sa:
                    # Check if the box is large enough to hold the ball cluster
                    # The ball radius is 2, so clearance is 4.
                    if l < 4 or w < 4 or h < 4:
                        continue

                    space_dims = sorted([l - 4, w - 4, h - 4])
                    
                    # Check if the cluster fits in the available space
                    if (space_dims[0] >= cluster_dims[0] and
                        space_dims[1] >= cluster_dims[1] and
                        space_dims[2] >= cluster_dims[2]):
                        
                        min_sa = sa
                        best_box = (l, w, h)

    if best_box:
        # Sort dimensions for consistent output, largest first
        l, w, h = sorted(best_box, reverse=True)
        d = int(min_sa)
        print(f"{l}:{w}:{h}:{d}")
    else:
        print("0")

solve()
