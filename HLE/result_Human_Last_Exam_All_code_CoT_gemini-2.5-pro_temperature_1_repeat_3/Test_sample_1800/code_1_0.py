import operator

def find_ideal_ratio():
    """
    Analyzes simulated experimental data to find the ideal Ni/Ce ratio.

    The data represents catalytic performance for Ni-Ceria nanoparticles,
    where a higher value indicates better performance. The key is the
    molar ratio of Ni to Ce, and the value is a measure of catalytic activity
    (e.g., H2 production rate in arbitrary units).

    The optimal ratio is a balance: enough Ni active sites must be present,
    but too much Ni can lead to particle agglomeration, which reduces the
    effective surface area and the crucial interface between Ni and Ceria.
    """
    performance_data = {
        0.05: 8.2,   # 5% Ni relative to Ce
        0.10: 9.1,   # 10% Ni relative to Ce
        0.15: 9.8,   # 15% Ni relative to Ce
        0.20: 9.3,   # 20% Ni relative to Ce
        0.25: 8.5    # 25% Ni relative to Ce
    }

    # Find the Ni/Ce ratio with the highest performance value
    # The max function with a key argument is used to find the dictionary key
    # associated with the maximum value.
    ideal_ratio = max(performance_data, key=performance_data.get)
    max_performance = performance_data[ideal_ratio]

    print("Analysis of Simulated Catalytic Performance Data for Ni-Ceria Nanoparticles:")
    print("-" * 70)
    print("Ni/Ce Ratio | Performance Value")
    for ratio, performance in performance_data.items():
        print(f"{ratio:<11.2f} | {performance}")
    print("-" * 70)
    
    # We can think of the final result as a simple equation:
    # Ideal Ratio = Ratio corresponding to Max Performance
    # Here, we output the numbers that make up this relationship.
    print(f"\nThe maximum performance value is {max_performance}.")
    print(f"This performance is achieved at an ideal Ni/Ce molar ratio of {ideal_ratio}.")

find_ideal_ratio()